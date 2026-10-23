#!/usr/bin/env python3
"""Export a Hail Variant Dataset as a filtered, annotated dense MatrixTable.

The VDS is densified with ``hl.vds.to_dense_mt``, so samples covered by a reference block at a variant site
get a hom-ref entry. Multi-allelic sites are kept (not split). Entry fields (GT, AD, PL, ...) are built from
the VDS local fields and site-level ``info`` fields are aggregated as in ``annotate_vds.py``.
"""

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
from typing import Literal
from urllib.parse import urlparse

import hail as hl

LOGGER = logging.getLogger("vds-to-mt")

INFO_PREFIX = "info."
LOCAL_FIELDS = ("LA", "LGT", "LAD", "LPL", "LPGT")
ALT_ALLELES_RANGE_FIELD = "alt_alleles_range_array"
# Temporary row fields holding the raw aggregation results before the final `info` struct is selected.
SITE_INFO_FIELD = "_site_info"
AS_INFO_FIELD = "_as_info"
CALL_STATS_FIELD = "_call_stats"
QD_FIELD = "_qd"
DB_FIELD = "_db"
# GATK FORMAT/SB: [ref_fwd, ref_rev, alt_fwd, alt_rev]. Some callers write a fwd/rev pair per allele instead.
GATK_SB_LENGTH = 4


@dataclass(slots=True, frozen=True)
class Annotation:
    name: str  # "AD" for an entry field, "info.QD" for a row info field
    kind: Literal["entry", "info"]
    sources: tuple[str, ...]  # required variant_data entry fields (top-level or inside gvcf_info)
    group: str  # builder that computes it

    @property
    def field(self) -> str:
        return self.name.removeprefix(INFO_PREFIX)


ANNOTATIONS: tuple[Annotation, ...] = (
    Annotation("GT", "entry", ("LGT", "LA"), "local_to_global"),
    Annotation("AD", "entry", ("LAD", "LA"), "local_to_global"),
    Annotation("PL", "entry", ("LPL", "LA"), "local_to_global"),
    Annotation("PGT", "entry", ("LPGT", "LA"), "local_to_global"),
    Annotation("DP", "entry", ("DP",), "passthrough"),
    Annotation("GQ", "entry", ("GQ",), "passthrough"),
    Annotation("MIN_DP", "entry", ("MIN_DP",), "passthrough"),
    Annotation("PID", "entry", ("PID",), "passthrough"),
    Annotation("PS", "entry", ("PS",), "passthrough"),
    Annotation("RGQ", "entry", ("RGQ",), "passthrough"),
    Annotation("SB", "entry", ("SB",), "passthrough"),
    Annotation("info.DP", "info", ("DP",), "site_info"),
    Annotation("info.ReadPosRankSum", "info", ("ReadPosRankSum",), "site_info"),
    Annotation("info.MQRankSum", "info", ("MQRankSum",), "site_info"),
    Annotation("info.BaseQRankSum", "info", ("BaseQRankSum",), "site_info"),
    Annotation("info.MQ", "info", ("RAW_MQandDP",), "site_info"),
    Annotation("info.FS", "info", ("SB",), "site_info"),
    Annotation("info.SOR", "info", ("SB",), "site_info"),
    Annotation("info.AS_MQ", "info", ("AS_RAW_MQ", "LAD", "LA"), "as_info"),
    Annotation("info.AS_ReadPosRankSum", "info", ("AS_RAW_ReadPosRankSum", "LA"), "as_info"),
    Annotation("info.AS_MQRankSum", "info", ("AS_RAW_MQRankSum", "LA"), "as_info"),
    Annotation("info.AS_QD", "info", ("AS_QD", "LA"), "as_info"),
    Annotation("info.AS_FS", "info", ("AS_SB_TABLE", "LA"), "as_info"),
    Annotation("info.AS_SOR", "info", ("AS_SB_TABLE", "LA"), "as_info"),
    Annotation("info.DB", "info", ("DB",), "db"),
    # QD needs global PL, which the qd builder computes itself from LPL/LA (see qd_expr).
    Annotation("info.QD", "info", ("LPL", "LA", "LGT", "DP"), "qd"),
    Annotation("info.AN", "info", ("LGT", "LA"), "call_stats"),
    Annotation("info.AC", "info", ("LGT", "LA"), "call_stats"),
    Annotation("info.AF", "info", ("LGT", "LA"), "call_stats"),
)
ANNOTATIONS_BY_NAME = {annotation.name: annotation for annotation in ANNOTATIONS}
DEFAULT_ANNOTATIONS = tuple(annotation.name for annotation in ANNOTATIONS)

# Inputs of the gnomAD helpers: annotation name -> (helper aggregation argument, input field).
# The helper output field is the annotation name without "info.".
SITE_INFO_INPUTS: dict[str, tuple[str, str]] = {
    "info.DP": ("int32_sum_agg_fields", "DP"),
    "info.ReadPosRankSum": ("median_agg_fields", "ReadPosRankSum"),
    "info.MQRankSum": ("median_agg_fields", "MQRankSum"),
    "info.BaseQRankSum": ("median_agg_fields", "BaseQRankSum"),
    "info.MQ": ("array_sum_agg_fields", "RAW_MQandDP"),  # gnomAD turns RAW_MQandDP into MQ
    "info.FS": ("array_sum_agg_fields", "SB"),  # gnomAD derives FS and SOR from SB
    "info.SOR": ("array_sum_agg_fields", "SB"),
}
AS_INFO_INPUTS: dict[str, tuple[str, str]] = {
    "info.AS_MQ": ("sum_agg_fields", "AS_RAW_MQ"),  # gnomAD turns AS_RAW_MQ into AS_MQ
    "info.AS_ReadPosRankSum": ("median_agg_fields", "AS_RAW_ReadPosRankSum"),
    "info.AS_MQRankSum": ("median_agg_fields", "AS_RAW_MQRankSum"),
    "info.AS_QD": ("median_agg_fields", "AS_QD"),
    "info.AS_FS": ("array_sum_agg_fields", "AS_SB_TABLE"),  # gnomAD derives AS_FS and AS_SOR from AS_SB_TABLE
    "info.AS_SOR": ("array_sum_agg_fields", "AS_SB_TABLE"),
}


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
    annotations: tuple[str, ...]
    strict_annotations: bool
    fill_ref_dp: bool


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def annotation_list(value: str) -> tuple[str, ...]:
    names = tuple(dict.fromkeys(name.strip() for name in value.split(",") if name.strip()))
    if not names:
        raise argparse.ArgumentTypeError("expected at least one annotation name")
    unknown = [name for name in names if name not in ANNOTATIONS_BY_NAME]
    if unknown:
        raise argparse.ArgumentTypeError(
            f"unknown annotation(s): {', '.join(unknown)}. Use --list-annotations to see the available names."
        )
    return names


def format_annotation_table() -> str:
    rows = [("NAME", "KIND", "BUILDER", "REQUIRED SOURCES")]
    rows += [(a.name, a.kind, a.group, ", ".join(a.sources)) for a in ANNOTATIONS]
    widths = [max(len(row[i]) for row in rows) for i in range(3)]
    return "\n".join(
        f"{name:<{widths[0]}}  {kind:<{widths[1]}}  {group:<{widths[2]}}  {sources}"
        for name, kind, group, sources in rows
    )


class ListAnnotationsAction(argparse.Action):
    """Print the annotation registry and exit, without requiring positional arguments or starting Hail."""

    def __call__(self, parser: argparse.ArgumentParser, *_args: object, **_kwargs: object) -> None:
        print(format_annotation_table())
        print(
            "\nSources are looked up as top-level variant_data entry fields or inside gvcf_info "
            "(local fields LA/LGT/LAD/LPL/LPGT must be top-level)."
        )
        parser.exit()


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
        LOGGER.info(
            f"=== Total variants in VDS after filtering to <= {max_alleles} alleles: "
            f"{filtered_vds.variant_data.count_rows()}"
        )
    return filtered_vds


def filter_vds_to_non_ref_sites(vds: hl.vds.VariantDataset, count_variants: bool) -> hl.vds.VariantDataset:
    # Densifying takes reference blocks from reference_data, so dropping variant rows without a non-ref LGT
    # before to_dense_mt gives the same result as filtering the dense MatrixTable, with less work.
    variant_data = vds.variant_data
    variant_data = variant_data.filter_rows(hl.agg.any(variant_data.LGT.is_non_ref()))
    if count_variants:
        LOGGER.info(f"=== Total variants in VDS after non-reference filter: {variant_data.count_rows()}")
    return hl.vds.VariantDataset(vds.reference_data, variant_data)


def missing_sources(variant_data: hl.MatrixTable, annotation: Annotation) -> list[str]:
    entry_fields = set(variant_data.entry)
    gvcf_info_fields = set(variant_data.gvcf_info) if "gvcf_info" in variant_data.entry else set()
    return [
        source
        for source in annotation.sources
        if source not in entry_fields and (source in LOCAL_FIELDS or source not in gvcf_info_fields)
    ]


def resolve_annotations(
    variant_data: hl.MatrixTable,
    requested: Sequence[str],
    strict: bool,
) -> tuple[Annotation, ...]:
    """Return the requested annotations whose sources exist, in registry order. GT is always included."""
    wanted = {"GT", *requested}
    resolved: list[Annotation] = []
    errors: list[str] = []
    for annotation in ANNOTATIONS:
        if annotation.name not in wanted:
            continue
        missing = missing_sources(variant_data, annotation)
        if not missing:
            resolved.append(annotation)
        elif strict or annotation.name == "GT":
            errors.append(f"{annotation.name} (missing: {', '.join(missing)})")
        else:
            LOGGER.warning("Skipping annotation %s: missing source field(s) %s", annotation.name, ", ".join(missing))
    if errors:
        raise ValueError(
            "Cannot build the requested annotation(s), source fields are not in the VDS variant data: "
            + "; ".join(errors)
        )
    return tuple(resolved)


def drop_non_gatk_sb_annotations(
    variant_data: hl.MatrixTable,
    annotations: tuple[Annotation, ...],
    strict: bool,
) -> tuple[Annotation, ...]:
    """Drop every SB-based annotation unless all SB entries in the data follow the GATK 4-value notation."""
    sb_annotations = [a for a in annotations if "SB" in a.sources]
    if not sb_annotations:
        return annotations
    names = ", ".join(a.name for a in sb_annotations)
    LOGGER.info("=== Checking that all SB entries follow GATK notation (%d values): %s", GATK_SB_LENGTH, names)
    sb = entry_source(variant_data, "SB")
    sb_lengths = variant_data.aggregate_entries(hl.agg.filter(hl.is_defined(sb), hl.agg.counter(hl.len(sb))))
    non_gatk = {length: count for length, count in sb_lengths.items() if length != GATK_SB_LENGTH}
    if not non_gatk:
        # All SB annotations follow GATK format
        return annotations
    detail = ", ".join(f"{count} with {length} values" for length, count in sorted(non_gatk.items()))
    message = f"SB entries not in GATK notation ({detail}), cannot build {names}"
    if strict:
        raise ValueError(message)
    LOGGER.warning("Skipping annotation(s): %s", message)
    return tuple(a for a in annotations if a not in sb_annotations)


def entry_source(mt: hl.MatrixTable, field: str) -> hl.Expression:
    # Same lookup rule as gnomAD's _get_info_agg_expr: top-level entry field wins over gvcf_info.
    if field in mt.entry:
        return mt[field]
    return mt.gvcf_info[field]


def local_to_global_expr(mt: hl.MatrixTable, name: str) -> hl.Expression:
    n_alleles = hl.len(mt.alleles)
    if name == "GT":
        return hl.vds.lgt_to_gt(mt.LGT, mt.LA)
    if name == "PGT":
        return hl.vds.lgt_to_gt(mt.LPGT, mt.LA)
    if name == "AD":
        return hl.vds.local_to_global(mt.LAD, mt.LA, n_alleles, fill_value=0, number="R")
    if name == "PL":
        return hl.vds.local_to_global(mt.LPL, mt.LA, n_alleles, fill_value=hl.missing(hl.tint32), number="G")
    raise ValueError(f"No local-to-global builder for {name}")


def annotate_genotype_entries(mt: hl.MatrixTable, annotations: Sequence[Annotation]) -> hl.MatrixTable:
    names = [a.name for a in annotations if a.group == "local_to_global"]
    LOGGER.info("=== Building entry annotations from local fields: %s", ", ".join(names))
    return mt.annotate_entries(**{name: local_to_global_expr(mt, name) for name in names})


def fill_ref_dp(mt: hl.MatrixTable) -> hl.MatrixTable:
    if "DP" not in mt.entry or "MIN_DP" not in mt.entry:
        LOGGER.warning("Skipping --fill-ref-dp: DP and MIN_DP entry fields are both required")
        return mt
    LOGGER.info("=== Filling missing DP with MIN_DP")
    return mt.annotate_entries(DP=hl.coalesce(mt.DP, mt.MIN_DP))


def gnomad_agg_fields(annotations: Sequence[Annotation], inputs: dict[str, tuple[str, str]]) -> dict[str, list[str]]:
    agg_fields: dict[str, list[str]] = {
        "sum_agg_fields": [],
        "int32_sum_agg_fields": [],
        "median_agg_fields": [],
        "array_sum_agg_fields": [],
    }
    for annotation in annotations:
        argument, field = inputs[annotation.name]
        if field not in agg_fields[argument]:
            agg_fields[argument].append(field)
    return agg_fields


def site_info_expr(mt: hl.MatrixTable, annotations: Sequence[Annotation]) -> hl.StructExpression:
    from gnomad.utils.sparse_mt import get_site_info_expr

    agg_fields = gnomad_agg_fields(annotations, SITE_INFO_INPUTS)
    if [a.name for a in annotations] == ["info.DP"]:
        # gnomAD wraps the non-DP fields in hl.agg.filter, which Hail rejects when that struct is empty.
        # Same as its int32 DP sum over all genotypes.
        return hl.struct(DP=hl.int32(hl.agg.sum(entry_source(mt, "DP"))))
    return get_site_info_expr(
        mt,
        sum_agg_fields=agg_fields["sum_agg_fields"],
        int32_sum_agg_fields=agg_fields["int32_sum_agg_fields"],
        median_agg_fields=agg_fields["median_agg_fields"],
        array_sum_agg_fields=agg_fields["array_sum_agg_fields"],
    )


def as_info_expr(mt: hl.MatrixTable, annotations: Sequence[Annotation]) -> hl.StructExpression:
    from gnomad.utils.sparse_mt import get_as_info_expr

    agg_fields = gnomad_agg_fields(annotations, AS_INFO_INPUTS)
    return get_as_info_expr(
        mt,
        sum_agg_fields=agg_fields["sum_agg_fields"],
        int32_sum_agg_fields=agg_fields["int32_sum_agg_fields"],
        median_agg_fields=agg_fields["median_agg_fields"],
        array_sum_agg_fields=agg_fields["array_sum_agg_fields"],
        alt_alleles_range_array_field=ALT_ALLELES_RANGE_FIELD,
        treat_fields_as_allele_specific=True,
    )


def qd_expr(mt: hl.MatrixTable) -> hl.Float64Expression:
    # Port of generate_qc_annotations from annotate_vds.py, aggregated in place on the dense rows.
    # PL is computed here from LPL/LA rather than taken from the PL entry annotation, so info.QD does not force
    # PL into the output: the PL entry field is only kept when it was requested.
    pl = hl.vds.local_to_global(mt.LPL, mt.LA, hl.len(mt.alleles), fill_value=hl.missing(hl.tint32), number="G")
    lin = 10 ** (-pl / 10.0)
    gp0 = lin[0] / hl.sum(lin)
    non_ref = mt.LGT.is_non_ref()
    qual = -10 * hl.agg.filter(non_ref, hl.agg.sum(hl.if_else(pl[0] > 3000, -300, hl.log10(gp0))))
    non_ref_dp = hl.agg.filter(non_ref, hl.agg.sum(entry_source(mt, "DP")))
    # Missing instead of inf/NaN when no depth is recorded on the non-ref genotypes.
    return hl.or_missing(non_ref_dp > 0, hl.min(35, qual / non_ref_dp))


def info_value(mt: hl.MatrixTable, annotation: Annotation) -> hl.Expression:
    field = annotation.field
    if annotation.group == "site_info":
        value = mt[SITE_INFO_FIELD][field]
    elif annotation.group == "as_info":
        value = mt[AS_INFO_FIELD][field]
    elif annotation.group == "call_stats":
        value = mt[CALL_STATS_FIELD][field]
        if field in ("AC", "AF"):
            value = value[1:]  # one value per alt allele
    elif annotation.group == "qd":
        value = mt[QD_FIELD]
    elif annotation.group == "db":
        value = mt[DB_FIELD]
    else:
        raise ValueError(f"No info builder for {annotation.name}")

    # Fisher's exact test gives NaN without strand-bias counts, store it as missing as annotate_vds.py does.
    if field == "FS":
        return hl.or_missing(~hl.is_nan(value), value)
    if field == "AS_FS":
        return value.map(lambda x: hl.or_missing(~hl.is_nan(x), x))
    return value


def annotate_info(mt: hl.MatrixTable, annotations: Sequence[Annotation]) -> hl.MatrixTable:
    info_annotations = [a for a in annotations if a.kind == "info"]
    if not info_annotations:
        return mt
    LOGGER.info("=== Building row info annotations: %s", ", ".join(a.name for a in info_annotations))

    def in_group(group: str) -> list[Annotation]:
        return [a for a in info_annotations if a.group == group]

    if in_group("as_info"):
        mt = mt.annotate_rows(**{ALT_ALLELES_RANGE_FIELD: hl.range(1, hl.len(mt.alleles))})

    # First aggregate each builder group once into a temporary row field, then pick the requested fields.
    # This also drops gnomAD intermediates such as SB and AS_SB_TABLE.
    aggregations: dict[str, hl.Expression] = {}
    if in_group("site_info"):
        aggregations[SITE_INFO_FIELD] = site_info_expr(mt, in_group("site_info"))
    if in_group("as_info"):
        aggregations[AS_INFO_FIELD] = as_info_expr(mt, in_group("as_info"))
    if in_group("call_stats"):
        aggregations[CALL_STATS_FIELD] = hl.agg.call_stats(mt.GT, mt.alleles)
    if in_group("qd"):
        aggregations[QD_FIELD] = qd_expr(mt)
    if in_group("db"):
        aggregations[DB_FIELD] = hl.agg.any(hl.coalesce(entry_source(mt, "DB"), False))
    mt = mt.annotate_rows(**aggregations)

    mt = mt.annotate_rows(info=hl.struct(**{a.field: info_value(mt, a) for a in info_annotations}))
    return mt.drop(*aggregations)


def promote_passthrough_entries(mt: hl.MatrixTable, annotations: Sequence[Annotation]) -> hl.MatrixTable:
    # Copy only the requested gvcf_info fields; a top-level entry field of the same name is kept as is.
    promoted = [a.name for a in annotations if a.group == "passthrough" and a.name not in mt.entry]
    if not promoted:
        return mt
    LOGGER.info("=== Promoting gvcf_info field(s) to entries: %s", ", ".join(promoted))
    return mt.annotate_entries(**{name: mt.gvcf_info[name] for name in promoted})


def drop_existing_fields(mt: hl.MatrixTable, fields: Sequence[str]) -> hl.MatrixTable:
    fields_to_drop = [field for field in fields if field in mt.row or field in mt.entry]
    if not fields_to_drop:
        return mt
    LOGGER.info("=== Dropping MatrixTable field(s): %s", ", ".join(fields_to_drop))
    return mt.drop(*fields_to_drop)


def matrix_table_from_vds(config: Config) -> hl.MatrixTable:
    LOGGER.info("=== Reading VDS: %s", config.input_vds)
    vds = hl.vds.read_vds(hail_path(config.input_vds))

    vds = filter_vds_by_allele_count(vds, config.max_alleles, config.count_variants)

    LOGGER.info("=== Filtering to sites with at least one non-reference genotype")
    vds = filter_vds_to_non_ref_sites(vds, config.count_variants)

    LOGGER.info("=== Checking source fields for annotations")
    annotations = resolve_annotations(vds.variant_data, config.annotations, config.strict_annotations)
    annotations = drop_non_gatk_sb_annotations(vds.variant_data, annotations, config.strict_annotations)

    LOGGER.info("=== Creating dense MatrixTable")
    mt = hl.vds.to_dense_mt(vds)
    if config.count_variants:
        LOGGER.info(f"=== Total variants in dense MT: {mt.count_rows()}")

    if config.fill_ref_dp:
        mt = fill_ref_dp(mt)

    mt = annotate_genotype_entries(mt, annotations)
    # Info aggregations need the local fields and gvcf_info, so they run before those are dropped.
    mt = annotate_info(mt, annotations)
    mt = promote_passthrough_entries(mt, annotations)

    # Keep only the built entry annotations; this drops the local fields (LA, LGT, LAD, LPL, LPGT).
    # gvcf_info is kept here and dropped below unless --keep-gvcf-info was given.
    entry_fields = [a.name for a in annotations if a.kind == "entry"]
    if "gvcf_info" in mt.entry:
        entry_fields.append("gvcf_info")
    mt = mt.select_entries(*entry_fields)
    mt = drop_existing_fields(mt, [ALT_ALLELES_RANGE_FIELD])
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
        description=(
            "Export a Hail VDS to a filtered, annotated dense Hail MatrixTable. "
            "Multi-allelic sites are kept (not split)."
        ),
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
        help="Optional reference FASTA path or URI to register with Hail. Not required by the export itself.",
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
    annotation_group = parser.add_mutually_exclusive_group()
    annotation_group.add_argument(
        "--annotations",
        type=annotation_list,
        default=None,
        metavar="FIELD[,FIELD...]",
        help=(
            "Comma-separated annotations to build. Entry fields are plain names (GT, AD), row INFO fields use "
            "the prefix 'info.' (info.QD). Annotations given here fail if their source fields are missing. "
            "Default: all annotations, skipping those whose source fields are missing."
        ),
    )
    annotation_group.add_argument(
        "--no-annotations",
        action="store_true",
        help="Build no annotations except GT, which is always built.",
    )
    parser.add_argument(
        "--strict-annotations",
        action="store_true",
        help="Fail if a default annotation cannot be built instead of logging a warning and skipping it.",
    )
    parser.add_argument(
        "--list-annotations",
        action=ListAnnotationsAction,
        nargs=0,
        default=argparse.SUPPRESS,
        help="Print the available annotations with their required source fields and exit.",
    )
    parser.add_argument(
        "--fill-ref-dp",
        action="store_true",
        help="Fill missing entry DP with MIN_DP (gnomAD convention). Changes info.DP. Default: False.",
    )
    parser.add_argument(
        "--keep-gvcf-info",
        action="store_true",
        help="Keep the gvcf_info entry field in the output MatrixTable. Can cause issues with exporting VCFs.",
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

    if args.no_annotations:
        annotations: tuple[str, ...] = ()
    elif args.annotations is not None:
        annotations = args.annotations
    else:
        annotations = DEFAULT_ANNOTATIONS

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
        annotations=annotations,
        # An explicitly requested list is always strict.
        strict_annotations=args.strict_annotations or args.annotations is not None,
        fill_ref_dp=args.fill_ref_dp,
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
