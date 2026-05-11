#!/usr/bin/env python3
"""Combine a folder of gVCFs into a Hail Variant Dataset."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import hail as hl

GVCF_SUFFIXES = (".g.vcf.gz", ".g.vcf.bgz")
LOGGER = logging.getLogger("vds-from-gvcf")

@dataclass(slots=True, frozen=True)
class Config:
    gvcf_dir: Path
    temp_vds_dir: str
    output_vds: str
    shard_size: int
    gvcf_batch_size: int
    reference: str
    tmp_dir: str | None
    whole_genome: bool
    recursive: bool
    overwrite: bool

def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def discover_gvcfs(gvcf_dir: Path, recursive: bool) -> list[str]:
    if not gvcf_dir.exists():
        raise FileNotFoundError(f"gVCF folder does not exist: {gvcf_dir}")
    if not gvcf_dir.is_dir():
        raise NotADirectoryError(f"gVCF path is not a folder: {gvcf_dir}")

    files: Iterable[Path]
    files = gvcf_dir.rglob("*") if recursive else gvcf_dir.iterdir()
    gvcfs = sorted(
        str(path.resolve())
        for path in files
        if path.is_file() and path.name.lower().endswith(GVCF_SUFFIXES)
    )
    if not gvcfs:
        suffixes = ", ".join(GVCF_SUFFIXES)
        raise ValueError(f"No gVCF files found in {gvcf_dir}. Expected suffixes: {suffixes}")
    return gvcfs


def chunked(items: Sequence[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(items), size):
        yield list(items[start : start + size])


def join_hail_path(base: str, *parts: str) -> str:
    normalized_base = base.rstrip("/")
    return "/".join([normalized_base, *(part.strip("/") for part in parts)])


def path_exists(path: str) -> bool:
    return Path(path).exists()


def remove_path(path: str) -> None:
    local_path = Path(path)
    if local_path.is_dir():
        shutil.rmtree(local_path)
    elif local_path.exists():
        local_path.unlink()

def prepare_output_path(path: str, overwrite: bool) -> None:
    if not path_exists(path):
        return
    if not overwrite:
        LOGGER.warning(f"Output path already exists: {path}. Use --overwrite to replace it.")
        return
    LOGGER.info("Removing existing output path: %s", path)
    remove_path(path)

def init_hail(config: Config):
    tmp_dir = config.tmp_dir or config.temp_vds_dir
    LOGGER.info(f"=== Initializing Hail with reference={config.reference} tmp_dir={tmp_dir}")
    hl.init(tmp_dir=tmp_dir)
    hl.default_reference(config.reference)

def combine_gvcfs(
    config: Config,
    output_path: str,
    temp_path: str,
    save_path: str,
    gvcf_paths: Sequence[str],

) -> None:
    prepare_output_path(output_path, config.overwrite)
    prepare_output_path(save_path, config.overwrite)

    start = time.monotonic()
    LOGGER.info(f"=== Combining {len(gvcf_paths)} gVCFs into {output_path}")
    combiner = hl.vds.new_combiner(
        output_path=output_path,
        gvcf_paths=list(gvcf_paths),
        temp_path=temp_path,
        save_path=save_path,
        use_genome_default_intervals=config.whole_genome,
        use_exome_default_intervals=(not config.whole_genome),
        gvcf_batch_size=config.gvcf_batch_size,
    )
    combiner.run()
    LOGGER.info(f"Finished {output_path} in {(time.monotonic() - start):.1f} seconds")


def create_shard_vds(gvcfs: Sequence[str], config: Config) -> list[str]:
    shard_paths: list[str] = []
    shards = list(chunked(gvcfs, config.shard_size))
    LOGGER.info(f"=== Split {len(gvcfs)} gVCFs into {len(shards)} shard(s) ===")

    for index, shard_gvcfs in enumerate(shards, start=1):
        shard_name = f"shard-{index:05d}.vds"
        shard_path = join_hail_path(config.temp_vds_dir, "shards", shard_name)
        save_path = join_hail_path(config.temp_vds_dir, "combiner-state", f"shard-{index:05d}")
        combine_gvcfs(    # The combiner reuses existing VDS, so we can rerun it.
            output_path=shard_path,
            temp_path=config.temp_vds_dir,
            save_path=save_path,
            gvcf_paths=shard_gvcfs,
            config=config,
        )
        shard_paths.append(shard_path)
    return shard_paths


def merge_shards(shard_paths: list[str], config: Config) -> None:
    save_path = join_hail_path(config.temp_vds_dir, "combiner-state", "final-merge")

    LOGGER.info(f"=== Merging {len(shard_paths)} VDS shards into {config.output_vds} ===")
    start = time.monotonic()
    combiner = hl.vds.new_combiner(
        output_path=config.output_vds,
        vds_paths=shard_paths,
        save_path=save_path,
        temp_path=config.temp_vds_dir,
        use_genome_default_intervals=config.whole_genome,
        use_exome_default_intervals=(not config.whole_genome),
    )
    combiner.run()
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
        "temp_vds_dir",
        help="Temporary folder where per-shard and merge-intermediate VDS outputs are written.",
    )
    parser.add_argument(
        "output_vds",
        help="Final VDS output path.",
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
        "--reference",
        default="GRCh38",
        help="Hail reference genome. Default: GRCh38.",
    )
    parser.add_argument(
        "--tmp-dir",
        help="Hail/Spark temporary directory. Defaults to temp_vds_dir.",
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
        temp_vds_dir=args.temp_vds_dir,
        output_vds=args.output_vds,
        shard_size=args.shard_size,
        gvcf_batch_size=args.gvcf_batch_size,
        reference=args.reference,
        tmp_dir=args.tmp_dir,
        whole_genome=args.whole_genome,
        recursive=args.recursive,
        overwrite=args.overwrite,
    )

def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()

    config = parse_args(argv)
    gvcfs = discover_gvcfs(config.gvcf_dir, config.recursive)
    LOGGER.info(f"=== Discovered {len(gvcfs)} gVCF file(s)" )
    init_hail(config)
    LOGGER.info("=== Hail initialized" )
    shard_paths = create_shard_vds(gvcfs, config)
    merge_shards(shard_paths, config)

    LOGGER.info("VDS build completed successfully: %s", config.output_vds)
    return 0


if __name__ == "__main__":
    sys.exit(main())
