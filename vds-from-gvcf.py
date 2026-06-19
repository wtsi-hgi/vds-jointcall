#!/usr/bin/env python3
"""Combine a folder of gVCFs into a Hail Variant Dataset."""

import argparse
import gzip
import json
import logging
import os
import re
import shutil
import sys
import time
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import hail as hl

GVCF_SUFFIXES = (".g.vcf.gz", ".g.vcf.bgz")
DEFAULT_CALL_FIELDS = []
LOGGER = logging.getLogger("vds-from-gvcf")
VDS_SUCCESS_MARKER_SUFFIX = ".success"
MANIFEST_SCHEMA_VERSION = 2
SHARD_MANIFEST_TYPE = "vds-jointcall-shard"
FINAL_MANIFEST_TYPE = "vds-jointcall-final"
SHARD_NAME_RE = re.compile(r"^shard-(\d{5})\.vds(?:\.success)?$")


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


@dataclass(slots=True, frozen=True)
class GvcfInput:
    filename: str
    path: str
    sample_id: str
    size: int
    mtime: float


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


def read_gvcf_sample_id(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                if len(columns) != 10:
                    raise ValueError(f"Expected exactly one sample column in {path}, found {max(len(columns) - 9, 0)}")
                sample_id = columns[9].strip()
                if not sample_id:
                    raise ValueError(f"Empty sample ID in gVCF header: {path}")
                return sample_id
    raise ValueError(f"Could not find #CHROM header with sample ID in gVCF: {path}")


def gvcf_input_from_path(path: Path) -> GvcfInput:
    stat = path.stat()
    return GvcfInput(
        filename=path.name,
        path=str(path),
        sample_id=read_gvcf_sample_id(path),
        size=stat.st_size,
        mtime=stat.st_mtime,
    )


def build_gvcf_inputs(gvcf_paths: Sequence[Path]) -> list[GvcfInput]:
    inputs = [gvcf_input_from_path(path) for path in gvcf_paths]
    validate_unique_current_inputs(inputs)
    return inputs


def find_duplicates(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def validate_unique_current_inputs(gvcfs: Sequence[GvcfInput]) -> None:
    duplicate_sample_ids = find_duplicates(gvcf.sample_id for gvcf in gvcfs)
    duplicate_names = find_duplicates(gvcf.filename for gvcf in gvcfs)
    if duplicate_sample_ids:
        ValueError(f"Ambiguous gVCF sample IDs detected: {', '.join(duplicate_sample_ids)}")
    if duplicate_names:
        raise ValueError(f"Ambiguous gVCF names detected: {', '.join(duplicate_names)}")


def chunked[T](items: Sequence[T], size: int) -> Iterable[list[T]]:
    for start in range(0, len(items), size):
        yield list(items[start : start + size])


def vds_success_marker_path(path: Path) -> Path:
    return path.parent / f"{path.name}{VDS_SUCCESS_MARKER_SUFFIX}"


def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Manifest is not a JSON object: {path}")
    return data


def combiner_config_manifest(config: Config) -> dict[str, Any]:
    return {
        "reference": config.reference,
        "whole_genome": config.whole_genome,
        "call_fields": config.call_fields,
        "shard_size": config.shard_size,
        "schema_version": MANIFEST_SCHEMA_VERSION,
    }


def validate_manifest_config(manifest: dict[str, Any], config: Config, manifest_path: Path) -> None:
    manifest_config = manifest.get("combiner_config")
    if not isinstance(manifest_config, dict):
        raise ValueError(f"Manifest missing combiner_config: {manifest_path}")

    expected = combiner_config_manifest(config)
    for key in ("reference", "whole_genome", "call_fields", "schema_version"):
        if manifest_config.get(key) != expected[key]:
            raise ValueError(
                f"Incompatible manifest {manifest_path}: {key}={manifest_config.get(key)!r}, "
                f"current={expected[key]!r}"
            )


def mark_vds_success(path: Path, manifest: dict[str, Any]) -> None:
    marker_path = vds_success_marker_path(path)
    write_json(marker_path, manifest)


def vds_shard_completed(path: Path) -> bool:
    return vds_success_marker_path(path).exists()


def shard_path_from_manifest_path(manifest_path: Path) -> Path:
    if not manifest_path.name.endswith(VDS_SUCCESS_MARKER_SUFFIX):
        raise ValueError(f"Not a success manifest path: {manifest_path}")
    return manifest_path.with_name(manifest_path.name[: -len(VDS_SUCCESS_MARKER_SUFFIX)])


def shard_index(path: Path) -> int | None:
    match = SHARD_NAME_RE.match(path.name)
    if match is None:
        return None
    return int(match.group(1))


def next_shard_index(shards_dir: Path) -> int:
    indexes: list[int] = []
    if shards_dir.exists():
        for path in shards_dir.iterdir():
            index = shard_index(path)
            if index is not None:
                indexes.append(index)
    return max(indexes, default=0) + 1


def new_shard_path(config: Config, index: int) -> Path:
    return config.temp_vds_dir / "shards" / f"shard-{index:05d}.vds"


def shard_save_path(config: Config, index: int) -> Path:
    return config.temp_vds_dir / "combiner-state" / f"shard-{index:05d}"


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
    gvcfs: Sequence[GvcfInput],
) -> None:
    start = time.monotonic()
    LOGGER.info(f"=== Combining {len(gvcfs)} gVCFs into {output_path}")
    combiner = hl.vds.new_combiner(
        output_path=str(output_path),
        gvcf_paths=[gvcf.path for gvcf in gvcfs],
        temp_path=str(temp_path),
        save_path=str(save_path),
        use_genome_default_intervals=config.whole_genome,
        use_exome_default_intervals=(not config.whole_genome),
        gvcf_batch_size=config.gvcf_batch_size,
        call_fields=config.call_fields,
    )
    combiner.run()
    mark_vds_success(output_path, shard_manifest(output_path, gvcfs, config))
    LOGGER.info(f"Finished {output_path} in {(time.monotonic() - start):.1f} seconds")


def shard_manifest(shard_path: Path, gvcfs: Sequence[GvcfInput], config: Config) -> dict[str, Any]:
    return {
        "manifest_type": SHARD_MANIFEST_TYPE,
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "created_at": datetime.now(UTC).isoformat(),
        "shard_path": str(shard_path),
        "combiner_config": combiner_config_manifest(config),
        "gvcfs": [asdict(gvcf) for gvcf in gvcfs],
    }


def gvcf_input_from_manifest(data: dict[str, Any], manifest_path: Path) -> GvcfInput:
    required = ("filename", "path", "sample_id", "size", "mtime")
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Manifest entry in {manifest_path} missing fields: {', '.join(missing)}")
    try:
        return GvcfInput(
            filename=str(data["filename"]),
            path=str(data["path"]),
            sample_id=str(data["sample_id"]),
            size=int(data["size"]),
            mtime=float(data["mtime"]),
        )
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid gVCF entry in manifest {manifest_path}: {data!r}") from error


def load_shard_manifest(manifest_path: Path, config: Config) -> tuple[Path, list[GvcfInput]]:
    shard_path = shard_path_from_manifest_path(manifest_path)
    if not shard_path.exists():
        raise ValueError(f"Shard manifest exists but VDS path is missing: {manifest_path} -> {shard_path}")

    manifest = read_json(manifest_path)
    if manifest.get("manifest_type") != SHARD_MANIFEST_TYPE:
        raise ValueError(f"Invalid shard manifest type in {manifest_path}: {manifest.get('manifest_type')!r}")
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ValueError(f"Unsupported shard manifest schema in {manifest_path}: {manifest.get('schema_version')!r}")
    validate_manifest_config(manifest, config, manifest_path)

    manifest_gvcfs = manifest.get("gvcfs")
    if not isinstance(manifest_gvcfs, list) or not manifest_gvcfs:
        raise ValueError(f"Shard manifest must include at least one gVCF: {manifest_path}")
    gvcfs = []
    for item in manifest_gvcfs:
        if not isinstance(item, dict):
            raise ValueError(f"Invalid gVCF entry in manifest {manifest_path}: {item!r}")
        gvcfs.append(gvcf_input_from_manifest(item, manifest_path))
    validate_unique_current_inputs(gvcfs)
    return shard_path, gvcfs


def load_completed_shards(config: Config) -> list[tuple[Path, list[GvcfInput]]]:
    shards_dir = config.temp_vds_dir / "shards"
    if not shards_dir.exists():
        return []

    completed_shards: list[tuple[Path, list[GvcfInput]]] = []
    for manifest_path in sorted(shards_dir.glob(f"*.vds{VDS_SUCCESS_MARKER_SUFFIX}")):
        completed_shards.append(load_shard_manifest(manifest_path, config))
    return completed_shards


def remove_incomplete_shards(config: Config) -> None:
    shards_dir = config.temp_vds_dir / "shards"
    if not shards_dir.exists():
        return

    for shard_path in sorted(shards_dir.glob("*.vds")):
        if vds_shard_completed(shard_path):
            continue
        index = shard_index(shard_path)
        LOGGER.warning(f"!!! Removing existing incomplete shard: {shard_path}")
        remove_path(shard_path)
        if index is not None:
            remove_path(shard_save_path(config, index))


def validate_completed_shards(
    current_gvcfs: Sequence[GvcfInput],
    completed_shards: Sequence[tuple[Path, list[GvcfInput]]],
) -> tuple[list[GvcfInput], list[GvcfInput]]:
    current_by_sample = {gvcf.sample_id: gvcf for gvcf in current_gvcfs}
    current_by_name = {gvcf.filename: gvcf for gvcf in current_gvcfs}
    completed_by_sample: dict[str, GvcfInput] = {}
    completed_by_name: dict[str, GvcfInput] = {}

    for _shard_path, shard_gvcfs in completed_shards:
        for gvcf in shard_gvcfs:
            if gvcf.sample_id in completed_by_sample:
                raise ValueError(f"Duplicate sample ID across completed shards: {gvcf.sample_id}")
            if gvcf.filename in completed_by_name:
                raise ValueError(f"Duplicate gVCF file name across completed shards: {gvcf.filename}")
            completed_by_sample[gvcf.sample_id] = gvcf
            completed_by_name[gvcf.filename] = gvcf

            current_by_same_sample = current_by_sample.get(gvcf.sample_id)
            current_by_same_name = current_by_name.get(gvcf.filename)
            if current_by_same_sample is not None and current_by_same_sample.filename != gvcf.filename:
                raise ValueError(
                    f"Ambiguous rename for sample {gvcf.sample_id}: completed name "
                    f"{gvcf.filename}, current name {current_by_same_sample.filename}"
                )
            if current_by_same_name is not None and current_by_same_name.sample_id != gvcf.sample_id:
                raise ValueError(
                    f"Ambiguous file name reuse for {gvcf.filename}: completed sample "
                    f"{gvcf.sample_id}, current sample {current_by_same_name.sample_id}"
                )
            if current_by_same_sample is not None and (
                current_by_same_sample.size != gvcf.size or current_by_same_sample.mtime != gvcf.mtime
            ):
                raise ValueError(f"Changed-in-place gVCF detected for sample {gvcf.sample_id} ({gvcf.filename})")

    missing_from_current = [
        gvcf for gvcf in completed_by_sample.values() if gvcf.sample_id not in current_by_sample
    ]
    new_gvcfs = [gvcf for gvcf in current_gvcfs if gvcf.sample_id not in completed_by_sample]
    return new_gvcfs, missing_from_current


def prepare_overwrite(config: Config) -> None:
    remove_path(config.temp_vds_dir / "shards")
    remove_path(config.temp_vds_dir / "combiner-state")
    remove_path(config.output_vds)
    remove_path(vds_success_marker_path(config.output_vds))


def create_shard_vds(gvcfs: list[GvcfInput], config: Config) -> list[Path]:
    if config.overwrite:
        LOGGER.info("=== Removing existing combiner state because --overwrite was requested ===")
        prepare_overwrite(config)

    LOGGER.info(f"=== Found {len(gvcfs)} samples to combine ===")
    remove_incomplete_shards(config)
    completed_shards = load_completed_shards(config)
    LOGGER.info(f"=== Loaded {len(completed_shards)} completed shard(s) from previous run(s) ===")
    new_gvcfs, missing_from_current = validate_completed_shards(gvcfs, completed_shards)
    LOGGER.info(f"=== New gVCFs to combine: {len(new_gvcfs)} ===")
    for missing in missing_from_current:
        LOGGER.warning(
            "!!! Completed shard sample is absent from current input but will be included in final merge: "
            "%s (%s)",
            missing.sample_id,
            missing.filename,
        )

    shard_paths = [shard_path for shard_path, _ in completed_shards]
    shards = list(chunked(new_gvcfs, config.shard_size))
    LOGGER.info("=== Reusing %s completed shard(s) ===", len(completed_shards))
    LOGGER.info(f"=== Split {len(new_gvcfs)} new gVCFs into {len(shards)} shard(s) ===")

    next_index = next_shard_index(config.temp_vds_dir / "shards")
    for shard_gvcfs in shards:
        shard_path = new_shard_path(config, next_index)
        save_path = shard_save_path(config, next_index)
        next_index += 1

        combine_gvcfs(  # The combiner reuses existing VDS, so we can rerun it.
            output_path=shard_path,
            temp_path=config.temp_vds_dir,
            save_path=save_path,
            gvcfs=shard_gvcfs,
            config=config,
        )
        shard_paths.append(shard_path)
        if config.verbose:
            vds_stats(shard_path)

    return shard_paths


def merge_shards(shard_paths: list[Path], config: Config) -> None:
    save_path = config.temp_vds_dir / "combiner-state" / "final-merge"
    if save_path.exists():
        remove_path(save_path)
    remove_path(config.output_vds)
    remove_path(vds_success_marker_path(config.output_vds))

    LOGGER.info(f"=== Merging {len(shard_paths)} VDS shards into {config.output_vds} ===")
    start = time.monotonic()
    combiner = hl.vds.new_combiner(
        output_path=str(config.output_vds),
        vds_paths=[str(path) for path in shard_paths],
        save_path=str(save_path),
        temp_path=str(config.temp_vds_dir),
        use_genome_default_intervals=config.whole_genome,
        use_exome_default_intervals=(not config.whole_genome),
        call_fields=config.call_fields,
    )
    combiner.run()
    mark_vds_success(config.output_vds, final_manifest(shard_paths, config))
    LOGGER.info(f"Finished {config.output_vds} in {(time.monotonic() - start):.1f} seconds")


def final_manifest(shard_paths: Sequence[Path], config: Config) -> dict[str, Any]:
    shards = []
    sample_ids = []
    for shard_path in shard_paths:
        manifest_path = vds_success_marker_path(shard_path)
        manifest = read_json(manifest_path)
        gvcfs = manifest.get("gvcfs")
        if not isinstance(gvcfs, list):
            raise ValueError(f"Shard manifest missing gVCF list: {manifest_path}")
        shards.append({"shard_path": str(shard_path), "gvcfs": gvcfs})
        for gvcf in gvcfs:
            if not isinstance(gvcf, dict) or "sample_id" not in gvcf:
                raise ValueError(f"Invalid gVCF entry in shard manifest: {manifest_path}")
            sample_ids.append(str(gvcf["sample_id"]))

    return {
        "manifest_type": FINAL_MANIFEST_TYPE,
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "created_at": datetime.now(UTC).isoformat(),
        "output_vds": str(config.output_vds),
        "combiner_config": combiner_config_manifest(config),
        "shard_paths": [str(path) for path in shard_paths],
        "samples": sorted(sample_ids),
        "shards": shards,
    }


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
        default=50,
        help="Number of gVCFs per shard VDS.",
    )
    parser.add_argument(
        "--gvcf-batch-size",
        type=positive_int,
        default=25,
        help="gVCF batch size passed to Hail's combiner. Whne <= --shard-size, Hail runs its own batching internally.",
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
        call_fields=args.call_fields if len(args.call_fields) > 0 else DEFAULT_CALL_FIELDS,
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
    gvcf_inputs = build_gvcf_inputs(gvcfs)
    init_hail(config)
    LOGGER.info("=== Hail initialized")
    shard_paths = create_shard_vds(gvcf_inputs, config)
    merge_shards(shard_paths, config)

    LOGGER.info("VDS build completed successfully: %s", config.output_vds)
    return 0


if __name__ == "__main__":
    sys.exit(main())
