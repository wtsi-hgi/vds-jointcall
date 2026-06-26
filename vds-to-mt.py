#!/usr/bin/env python3
"""Export a Hail Variant Dataset as a filtered sparse MatrixTable."""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import hail as hl

LOGGER = logging.getLogger("vds-to-mt")


@dataclass(slots=True, frozen=True)
class Config:
    input_vds: str
    output_mt: str
    reference: str
    reference_sequence: str | None
    tmp_dir: str
    spark_memory: str
    max_alleles: int
    drop_gvcf_info: bool
    count_variants: bool
    overwrite: bool
    describe: bool


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def has_uri_scheme(path: str) -> bool:
    return urlparse(path).scheme != ""


def hail_path(path: str) -> str:
    if has_uri_scheme(path):
        return path
    return Path(path).resolve().as_uri()


def local_path(path: str) -> Path | None:
    if has_uri_scheme(path):
        return None
    return Path(path)


def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def init_hail(config: Config) -> None:
    if config.spark_memory != "":
        os.environ["PYSPARK_SUBMIT_ARGS"] = (
            f"--driver-memory {config.spark_memory} --executor-memory {config.spark_memory} pyspark-shell"
        )
        LOGGER.info("=== PYSPARK_SUBMIT_ARGS=%s", os.environ["PYSPARK_SUBMIT_ARGS"])
    else:
        LOGGER.info("=== Skipping PYSPARK_SUBMIT_ARGS because Spark memory option is empty")

    LOGGER.info("=== Initializing Hail with reference=%s tmp_dir=%s", config.reference, config.tmp_dir)
    hl.init(
        tmp_dir=hail_path(config.tmp_dir),
        spark_conf={
            'spark.rpc.message.maxSize': '1024',   # Required for VDS export for large datasets
            'spark.driver.maxResultSize': '0'  # 0 means unlimited; prevents driver OOMs during large collections
        }
    )
    hl.default_reference(config.reference)

    if config.reference_sequence is not None:
        LOGGER.info("=== Adding reference sequence for %s: %s", config.reference, config.reference_sequence)
        hl.get_reference(config.reference).add_sequence(hail_path(config.reference_sequence))


def validate_paths(config: Config) -> None:
    input_vds = local_path(config.input_vds)
    if input_vds is not None:
        if not input_vds.exists():
            raise FileNotFoundError(f"Input VDS does not exist: {input_vds}")
        if not input_vds.is_dir():
            raise NotADirectoryError(f"Input VDS path is not a directory: {input_vds}")

    output_mt = local_path(config.output_mt)
    if output_mt is not None:
        if output_mt.exists():
            if not config.overwrite:
                raise FileExistsError(f"Output MatrixTable already exists, use --overwrite to replace it: {output_mt}")
            LOGGER.info("=== Removing existing output because --overwrite was requested: %s", output_mt)
            remove_path(output_mt)
        output_mt.parent.mkdir(parents=True, exist_ok=True)

    tmp_dir = local_path(config.tmp_dir)
    if tmp_dir is not None:
        tmp_dir.mkdir(parents=True, exist_ok=True)


def filter_vds_by_allele_count(
    vds: hl.vds.VariantDataset,
    max_alleles: int,
    count_variants: bool,
) -> hl.vds.VariantDataset:
    variant_rows = vds.variant_data.rows()
    if count_variants:
        LOGGER.info(f"=== Total variants in VDS: {variant_rows.count()}")

    sites_to_keep = variant_rows.filter(hl.len(variant_rows.alleles) <= max_alleles)
    filtered_vds = hl.vds.filter_variants(vds, sites_to_keep)
    if count_variants:
        LOGGER.info(f"=== Total variants in VDS after filtering to <= {max_alleles} alleles: {filtered_vds.variant_data.count()}" )
    return filtered_vds

def annotate_variant_entries(vds: hl.vds.VariantDataset) -> hl.MatrixTable:
    variant_data = vds.variant_data.annotate_entries(
        AD=hl.vds.local_to_global(
            vds.variant_data.LAD,
            vds.variant_data.LA,
            n_alleles=hl.len(vds.variant_data.alleles),
            fill_value=0,
            number="R",
        )
    )
    variant_data = variant_data.annotate_entries(GT=hl.vds.lgt_to_gt(variant_data.LGT, variant_data.LA))
    return variant_data


def drop_existing_fields(mt: hl.MatrixTable, fields: Sequence[str]) -> hl.MatrixTable:
    fields_to_drop = [field for field in fields if field in mt.row]
    if not fields_to_drop:
        return mt
    LOGGER.info("=== Dropping MatrixTable row field(s): %s", ", ".join(fields_to_drop))
    return mt.drop(*fields_to_drop)


def matrix_table_from_vds(config: Config) -> hl.MatrixTable:
    LOGGER.info("=== Reading VDS: %s", config.input_vds)
    vds = hl.vds.read_vds(hail_path(config.input_vds))

    vds = filter_vds_by_allele_count(vds, config.max_alleles, config.count_variants)

    LOGGER.info("=== Annotating variant entries with AD and GT")
    variant_data = annotate_variant_entries(vds)

    LOGGER.info("=== Creating merged sparse MatrixTable")
    mt = hl.vds.to_merged_sparse_mt(hl.vds.VariantDataset(vds.reference_data, variant_data))

    LOGGER.info("=== Splitting multi-allelic sites")
    mt = hl.split_multi_hts(mt)
    if config.count_variants:
        LOGGER.info(f"=== Total variants in MT after split_multi_hts {mt.count_rows()}")

    LOGGER.info("=== Filtering to sites with at least one non-reference genotype")
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    if config.count_variants:
        LOGGER.info(f"=== Total variants in MT after non-reference filter {mt.count_rows()}")

    if config.drop_gvcf_info:
        mt = drop_existing_fields(mt, ["gvcf_info"])

    return mt


def write_matrix_table(config: Config) -> None:
    start = time.monotonic()
    mt = matrix_table_from_vds(config)

    if config.describe:
        mt.describe()

    LOGGER.info("=== Writing MatrixTable: %s", config.output_mt)
    mt.write(hail_path(config.output_mt), overwrite=config.overwrite)
    LOGGER.info("Finished %s in %.1f seconds", config.output_mt, time.monotonic() - start)


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_args(argv: Sequence[str] | None = None) -> Config:
    parser = argparse.ArgumentParser(
        description="Export a Hail VDS to a filtered sparse Hail MatrixTable.",
    )
    parser.add_argument(
        "input_vds",
        help="Input Hail VDS path or URI.",
    )
    parser.add_argument(
        "output_mt",
        help="Output Hail MatrixTable path or URI.",
    )
    parser.add_argument(
        "--reference",
        default="GRCh38",
        help="Hail reference genome. Default: GRCh38.",
    )
    parser.add_argument(
        "--reference-sequence",
        default=None,
        help="Additional reference FASTA path or URI to register with Hail",
    )
    parser.add_argument(
        "--tmp-dir",
        default="hail-tmp",
        help="Hail/Spark temporary directory path or URI. Default: hail-tmp.",
    )
    parser.add_argument(
        "-M",
        "--spark-memory",
        type=str,
        default="",
        help=(
            "Spark driver and executor memory used in PYSPARK_SUBMIT_ARGS. "
            "Use an empty value to skip setting it."
        ),
    )
    parser.add_argument(
        "--max-alleles",
        type=positive_int,
        default=10,
        help="Keep only variants with this many total alleles or fewer. Default: 10.",
    )
    parser.add_argument(
        "--keep-gvcf-info",
        action="store_true",
        help="Keep the gvcf_info row field in the output MatrixTable. Can cause issues with exporting VCFs. Default: False.",
    )
    parser.add_argument(
        "--count_variants",
        action="store_true",
        help="Log the number of variations after each step. Increases processing time on large datasets.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output MatrixTable.",
    )
    parser.add_argument(
        "--describe",
        action="store_true",
        help="Print the final MatrixTable schema before writing.",
    )
    args = parser.parse_args(argv)

    return Config(
        input_vds=args.input_vds,
        output_mt=args.output_mt,
        reference=args.reference,
        reference_sequence=args.reference_sequence,
        tmp_dir=args.tmp_dir,
        spark_memory=args.spark_memory,
        max_alleles=args.max_alleles,
        drop_gvcf_info=not args.keep_gvcf_info,
        count_variants=args.count_variants,
        overwrite=args.overwrite,
        describe=args.describe,
    )


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    config = parse_args(argv)
    validate_paths(config)
    init_hail(config)
    write_matrix_table(config)
    LOGGER.info("MatrixTable export completed successfully: %s", config.output_mt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
