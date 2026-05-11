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
from urllib.parse import urlparse

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
    use_exome_default_intervals: bool
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


def is_uri(path: str) -> bool:
    return bool(urlparse(path).scheme)


def path_exists(path: str) -> bool:
    if is_uri(path):
        return False
    return Path(path).exists()


def remove_path(path: str) -> None:
    if is_uri(path):
        raise ValueError(f"Cannot delete URI path without Hail filesystem support: {path}")

    local_path = Path(path)
    if local_path.is_dir():
        shutil.rmtree(local_path)
    elif local_path.exists():
        local_path.unlink()


def prepare_output_path(path: str, overwrite: bool) -> None:
    if not path_exists(path):
        return
    if not overwrite:
        raise FileExistsError(f"Output path already exists: {path}. Use --overwrite to replace it.")
    LOGGER.info("Removing existing output path: %s", path)
    remove_path(path)


def init_hail(config: Config):
    tmp_dir = config.tmp_dir or config.temp_vds_dir
    LOGGER.info("Initializing Hail with reference=%s tmp_dir=%s", config.reference, tmp_dir)
    hl.init(tmp_dir=tmp_dir)
    hl.default_reference(config.reference)

def run_combiner(
    config: Config,
    output_path: str,
    temp_path: str,
    save_path: str,
    gvcf_paths: Sequence[str] | None = None,
    vds_paths: Sequence[str] | None = None,

) -> None:
    prepare_output_path(output_path, config.overwrite)
    prepare_output_path(save_path, config.overwrite)

    start = time.monotonic()
    if gvcf_paths is not None:
        LOGGER.info("Combining %d gVCFs into %s", len(gvcf_paths), output_path)
        combiner = hl.vds.new_combiner(
            output_path=output_path,
            gvcf_paths=list(gvcf_paths),
            temp_path=temp_path,
            save_path=save_path,
            use_exome_default_intervals=config.use_exome_default_intervals,
            gvcf_batch_size=config.gvcf_batch_size,
        )
    elif vds_paths is not None:
        LOGGER.info("Merging %d VDS shards into %s", len(vds_paths), output_path)
        combiner = hl.vds.new_combiner(
            output_path=output_path,
            vds_paths=list(vds_paths),
            temp_path=temp_path,
            save_path=save_path,
            use_exome_default_intervals=config.use_exome_default_intervals,
        )
    else:
        raise ValueError("Either gvcf_paths or vds_paths must be provided")

    combiner.run()
    LOGGER.info("Finished %s in %.1f seconds", output_path, time.monotonic() - start)


def create_shard_vds(gvcfs: Sequence[str], config: Config) -> list[str]:
    shard_paths: list[str] = []
    shards = list(chunked(gvcfs, config.shard_size))
    LOGGER.info("Split %d gVCFs into %d shard(s)", len(gvcfs), len(shards))

    for index, shard_gvcfs in enumerate(shards, start=1):
        shard_name = f"shard-{index:05d}.vds"
        shard_path = join_hail_path(config.temp_vds_dir, "shards", shard_name)
        save_path = join_hail_path(config.temp_vds_dir, "combiner-state", f"shard-{index:05d}")
        run_combiner(
            output_path=shard_path,
            temp_path=config.temp_vds_dir,
            save_path=save_path,
            gvcf_paths=shard_gvcfs,
            config=config,
        )
        shard_paths.append(shard_path)

    return shard_paths


def merge_shards(hl, shard_paths: Sequence[str], config: Config) -> None:
    save_path = join_hail_path(config.temp_vds_dir, "combiner-state", "final-merge")
    run_combiner(
        hl,
        output_path=config.output_vds,
        temp_path=config.temp_vds_dir,
        save_path=save_path,
        vds_paths=shard_paths,
        config=config,
    )


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
        use_exome_default_intervals=not args.whole_genome,
        recursive=args.recursive,
        overwrite=args.overwrite,
    )

def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()

    config = parse_args(argv)
    gvcfs = discover_gvcfs(config.gvcf_dir, config.recursive)
    LOGGER.info("Discovered %d gVCF file(s)", len(gvcfs))
    init_hail(config)
    shard_paths = create_shard_vds(gvcfs, config)
    merge_shards(shard_paths, config)

    LOGGER.info("VDS build completed successfully: %s", config.output_vds)
    return 0


if __name__ == "__main__":
    sys.exit(main())
