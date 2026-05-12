#!/usr/bin/env python3
"""Combine a folder of gVCFs into a Hail Variant Dataset."""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import hail as hl

GVCF_SUFFIXES = (".g.vcf.gz", ".g.vcf.bgz")
LOGGER = logging.getLogger("vds-from-gvcf")
VDS_SUCCESS_MARKER_SUFFIX = ".success"


@dataclass(slots=True, frozen=True)
class Config:
    gvcf_dir: Path
    temp_vds_dir: Path
    output_vds: Path
    shard_size: int
    gvcf_batch_size: int
    call_fields: list[str]
    reference: str
    tmp_dir: Path | None
    spark_memory: str
    whole_genome: bool
    recursive: bool
    overwrite: bool
    verbose: bool = False


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def comma_separated_fields(value: str) -> list[str]:
    if value == "":
        return []
    fields = [field.strip() for field in value.split(":") if field.strip()]
    if not fields:
        raise argparse.ArgumentTypeError("must include at least one field")
    return fields


def memory_size(value: str) -> str:
    if value == "":
        raise argparse.ArgumentTypeError("must not be empty")
    return value


def discover_gvcfs(gvcf_dir: Path, recursive: bool) -> list[Path]:
    if not gvcf_dir.exists():
        raise FileNotFoundError(f"gVCF folder does not exist: {gvcf_dir}")
    if not gvcf_dir.is_dir():
        raise NotADirectoryError(f"gVCF path is not a folder: {gvcf_dir}")

    files: Iterable[Path]
    files = gvcf_dir.rglob("*") if recursive else gvcf_dir.iterdir()
    gvcfs = sorted(path.resolve() for path in files if path.is_file() and path.name.lower().endswith(GVCF_SUFFIXES))
    if not gvcfs:
        suffixes = ", ".join(GVCF_SUFFIXES)
        raise ValueError(f"No gVCF files found in {gvcf_dir}. Expected suffixes: {suffixes}")
    return gvcfs


def chunked[T](items: Sequence[T], size: int) -> Iterable[list[T]]:
    for start in range(0, len(items), size):
        yield list(items[start : start + size])


def combiner_call_fields_kwargs(config: Config) -> dict[str, list[str]]:
    if len(config.call_fields) == 0:
        return {}
    return {"call_fields": config.call_fields}


def vds_success_marker_path(path: Path) -> Path:
    return path.parent / f"{path.name}{VDS_SUCCESS_MARKER_SUFFIX}"


def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def mark_vds_success(path: Path) -> None:
    marker_path = vds_success_marker_path(path)
    marker_path.touch()


def vds_shard_completed(path: Path) -> bool:
    return vds_success_marker_path(path).exists()


def vds_stats(shard_path: Path) -> None:
    vds = hl.vds.read_vds(str(shard_path))
    print(vds.variant_data.count())  # Count variants (rows x samples)
    vds.variant_data.rows().describe()  # Schema of variant table
    vds.variant_data.entry.describe()


def init_hail(config: Config):
    tmp_dir = config.tmp_dir or config.temp_vds_dir
    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {config.spark_memory} --executor-memory {config.spark_memory} pyspark-shell"
    )
    LOGGER.info(f"=== Initializing Hail with reference={config.reference} tmp_dir={tmp_dir}")
    LOGGER.info(f"=== PYSPARK_SUBMIT_ARGS={os.environ['PYSPARK_SUBMIT_ARGS']}")
    hl.init(tmp_dir=str(tmp_dir))
    hl.default_reference(config.reference)


def combine_gvcfs(
    config: Config,
    output_path: Path,
    temp_path: Path,
    save_path: Path,
    gvcf_paths: Sequence[Path],
) -> None:
    start = time.monotonic()
    LOGGER.info(f"=== Combining {len(gvcf_paths)} gVCFs into {output_path}")
    combiner = hl.vds.new_combiner(
        output_path=str(output_path),
        gvcf_paths=[str(path) for path in gvcf_paths],
        temp_path=str(temp_path),
        save_path=str(save_path),
        use_genome_default_intervals=config.whole_genome,
        use_exome_default_intervals=(not config.whole_genome),
        gvcf_batch_size=config.gvcf_batch_size,
        **combiner_call_fields_kwargs(config),
    )
    combiner.run()
    mark_vds_success(output_path)
    LOGGER.info(f"Finished {output_path} in {(time.monotonic() - start):.1f} seconds")


def create_shard_vds(gvcfs: list[Path], config: Config) -> list[Path]:
    shard_paths: list[Path] = []
    shards = list(chunked(gvcfs, config.shard_size))
    LOGGER.info(f"=== Split {len(gvcfs)} gVCFs into {len(shards)} shard(s) ===")

    for index, shard_gvcfs in enumerate(shards, start=1):
        shard_name = f"shard-{index:05d}.vds"
        shard_path = config.temp_vds_dir / "shards" / shard_name
        save_path = config.temp_vds_dir / "combiner-state" / f"shard-{index:05d}"

        if shard_path.exists():
            if vds_shard_completed(shard_path):
                if config.overwrite:
                    LOGGER.info(f"=== Overwriting existing completed shard: {shard_path} ===")
                    remove_path(shard_path)
                    remove_path(save_path)
                else:
                    # The VDS and the corresponding success lock exist, so we can skip this shard.
                    LOGGER.info(f"=== Skipping existing completed shard: {shard_path} ===")
                    LOGGER.info("Use --overwrite to replace it.")
                    shard_paths.append(shard_path)
                    continue
            else:
                LOGGER.warning(f"!!! Removing existing incomplete shard: {shard_path}")
                remove_path(shard_path)
                remove_path(save_path)

        combine_gvcfs(  # The combiner reuses existing VDS, so we can rerun it.
            output_path=shard_path,
            temp_path=config.temp_vds_dir,
            save_path=save_path,
            gvcf_paths=shard_gvcfs,
            config=config,
        )
        shard_paths.append(shard_path)
        if config.verbose:
            vds_stats(shard_path)

    return shard_paths


def merge_shards(shard_paths: list[Path], config: Config) -> None:
    save_path = config.temp_vds_dir / "combiner-state" / "final-merge"
    # Final VDS combination always works with overwriting
    if save_path.exists():
        remove_path(save_path)

    LOGGER.info(f"=== Merging {len(shard_paths)} VDS shards into {config.output_vds} ===")
    start = time.monotonic()
    combiner = hl.vds.new_combiner(
        output_path=str(config.output_vds),
        vds_paths=[str(path) for path in shard_paths],
        save_path=str(save_path),
        temp_path=str(config.temp_vds_dir),
        use_genome_default_intervals=config.whole_genome,
        use_exome_default_intervals=(not config.whole_genome),
        **combiner_call_fields_kwargs(config),
    )
    combiner.run()
    mark_vds_success(config.output_vds)
    LOGGER.info(f"Finished {config.output_vds} in {(time.monotonic() - start):.1f} seconds")


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_args(argv: Sequence[str] | None = None) -> Config:
    parser = argparse.ArgumentParser(
        description="Combine all gVCFs in a local folder into a Hail VDS.",
    )
    parser.add_argument(
        "gvcf_dir",
        type=Path,
        help="Local folder containing gVCF files.",
    )
    parser.add_argument(
        "output_vds",
        type=Path,
        help="Final VDS output path.",
    )
    parser.add_argument(
        "--temp-vds-dir",
        type=Path,
        default=None,
        help="Temporary folder where per-shard and merge-intermediate VDS outputs are written.",
    )
    parser.add_argument(
        "--shard-size",
        type=positive_int,
        default=5,
        help="Number of gVCFs per shard VDS.",
    )
    parser.add_argument(
        "--gvcf-batch-size",
        type=positive_int,
        default=5,
        help="gVCF batch size passed to Hail's combiner.",
    )
    parser.add_argument(
        "--call-fields",
        type=comma_separated_fields,
        default="",  #  "GT:AVG_GQ:GQ:MIN_DP:MIN_GQ:PL",
        help="Colon-separated FORMAT call fields to include, passed to Hail's combiner call_fields option.",
    )
    parser.add_argument(
        "--reference",
        default="GRCh38",
        help="Hail reference genome. Default: GRCh38.",
    )
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        help="Hail/Spark temporary directory. Defaults to temp_vds_dir.",
    )
    parser.add_argument(
        "--spark-memory",
        type=memory_size,
        default="12G",
        help="Spark driver and executor memory used in PYSPARK_SUBMIT_ARGS. Default: 12G.",
    )
    parser.add_argument(
        "--whole-genome",
        action="store_true",
        help="Do not use Hail's exome default intervals.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search gvcf_dir recursively.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete existing shard, intermediate, and final VDS paths before writing.",
    )
    args = parser.parse_args(argv)

    return Config(
        gvcf_dir=args.gvcf_dir,
        temp_vds_dir=args.temp_vds_dir if args.temp_vds_dir is not None
            else args.gvcf_dir.with_suffix('.vds-combine'),
        output_vds=args.output_vds,
        shard_size=args.shard_size,
        gvcf_batch_size=args.gvcf_batch_size,
        call_fields=args.call_fields,
        reference=args.reference,
        tmp_dir=args.tmp_dir,
        spark_memory=args.spark_memory,
        whole_genome=args.whole_genome,
        recursive=args.recursive,
        overwrite=args.overwrite,
    )


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()

    config = parse_args(argv)
    gvcfs = discover_gvcfs(config.gvcf_dir, config.recursive)
    LOGGER.info(f"=== Discovered {len(gvcfs)} gVCF file(s)")
    init_hail(config)
    LOGGER.info("=== Hail initialized")
    shard_paths = create_shard_vds(gvcfs, config)
    merge_shards(shard_paths, config)

    LOGGER.info("VDS build completed successfully: %s", config.output_vds)
    return 0


if __name__ == "__main__":
    sys.exit(main())
