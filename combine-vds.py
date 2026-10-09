#!/usr/bin/env python3
"""Combine two existing Hail VDS (VariantDataset) files into one merged VDS.

Usage:
    combine_vds.py input1.vds input2.vds output.vds
"""

import argparse
import json
import logging
import os
import shutil
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Sequence

import hail as hl

LOGGER = logging.getLogger("combine-vds")

SCRIPT_VERSION = "2026-08-27.2-type-and-order-alignment"

VDS_SUCCESS_MARKER_SUFFIX = ".success"
MANIFEST_SCHEMA_VERSION = 2
FINAL_MANIFEST_TYPE = "vds-combine-final"


def configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def comma_separated_fields(value: str) -> list[str]:
    if value == "":
        return []
    fields = [field.strip() for field in value.split(":") if field.strip()]
    if not fields:
        raise argparse.ArgumentTypeError("must include at least one field")
    return fields


def hail_file_uri(path: Path) -> str:
    return Path(path).resolve().as_uri()


def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def vds_success_marker_path(path: Path) -> Path:
    return path.parent / f"{path.name}{VDS_SUCCESS_MARKER_SUFFIX}"


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


def combiner_config_manifest(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "reference": args.reference,
        "whole_genome": args.whole_genome,
        "call_fields": args.call_fields,
        "schema_version": MANIFEST_SCHEMA_VERSION,
    }


def validate_input_manifest(path: Path, args: argparse.Namespace, label: str) -> None:
    """If an input VDS has a .success manifest (e.g. from the gVCF combiner script),
    check its combiner config matches what we're about to request. This is what would
    have caught the earlier 'requested call fields: [PGT]' mismatch before running."""
    manifest_path = vds_success_marker_path(path)
    if not manifest_path.exists():
        LOGGER.info("=== No .success manifest found for %s (%s); skipping consistency check", label, path)
        return

    manifest = read_json(manifest_path)
    manifest_config = manifest.get("combiner_config")
    if not isinstance(manifest_config, dict):
        LOGGER.warning("!!! %s manifest has no combiner_config; skipping consistency check: %s", label, manifest_path)
        return

    expected = combiner_config_manifest(args)
    mismatches = {
        key: (manifest_config.get(key), expected[key])
        for key in ("reference", "whole_genome", "call_fields")
        if manifest_config.get(key) != expected[key]
    }
    if mismatches:
        details = ", ".join(f"{key}: existing={old!r} requested={new!r}" for key, (old, new) in mismatches.items())
        raise ValueError(f"{label} ({path}) was built with different combiner settings than requested: {details}")


def read_vds_sample_ids(path: Path) -> list[str]:
    vds = hl.vds.read_vds(hail_file_uri(path))
    return sorted(vds.variant_data.s.collect())


def mark_vds_success(path: Path, manifest: dict[str, Any]) -> None:
    write_json(vds_success_marker_path(path), manifest)


def final_manifest(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "manifest_type": FINAL_MANIFEST_TYPE,
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "created_at": datetime.now(UTC).isoformat(),
        "output_vds": str(args.output_vds),
        "combiner_config": combiner_config_manifest(args),
        "input_vds": [str(args.input1.resolve()), str(args.input2.resolve())],
        "samples": read_vds_sample_ids(args.output_vds),
    }


def vds_entry_field_order(mt: "hl.MatrixTable") -> list[str]:
    return list(mt.entry.dtype.fields)


def entry_field_types(mt: "hl.MatrixTable") -> dict[str, "hl.HailType"]:
    dtype = mt.entry.dtype
    return dict(zip(dtype.fields, dtype.types))


def resolve_canonical_type(field: str, type1: "hl.HailType", type2: "hl.HailType") -> "hl.HailType":
    if type1 == type2:
        return type1
    if {type1, type2} == {hl.tcall, hl.tstr}:
        # One input parsed this field (e.g. PGT) into a Call, the other left it as a raw
        # string - this happens when the two VDS were built with different --call-fields
        # settings during gVCF combining. Prefer the parsed `call` type and parse the
        # string side to match.
        return hl.tcall
    raise ValueError(
        f"Entry field {field!r} has incompatible types between input1 and input2 that can't be "
        f"safely reconciled: {type1} vs {type2}. This usually means the two VDS were built with "
        "different combiner settings (e.g. --call-fields) at gVCF-combine time. Resolve manually "
        "before combining."
    )


def align_matrix_table_schema(
    mt: "hl.MatrixTable", canonical_order: list[str], canonical_types: dict[str, "hl.HailType"]
) -> "hl.MatrixTable":
    current_types = entry_field_types(mt)
    cast_exprs = {}
    for field, target_type in canonical_types.items():
        current_type = current_types[field]
        if current_type == target_type:
            continue
        if target_type == hl.tcall and current_type == hl.tstr:
            cast_exprs[field] = hl.parse_call(mt[field])
        elif target_type == hl.tstr and current_type == hl.tcall:
            cast_exprs[field] = hl.str(mt[field])
        else:
            raise ValueError(f"Cannot safely cast entry field {field!r} from {current_type} to {target_type}")
    if cast_exprs:
        mt = mt.annotate_entries(**cast_exprs)
    return mt.select_entries(*canonical_order)


def align_vds_for_combiner(args: argparse.Namespace) -> tuple[str, str]:
    """Hail's combiner fails with 'must have the same row type' when two VDS have
    entry fields in a different order, or the same field with a different type
    (e.g. LPGT as `call` in one VDS and `str` in the other, from differing
    --call-fields settings at gVCF-combine time). Read both VDS and, if needed,
    rewrite copies with a canonical field order and reconciled types to
    temp_dir/aligned/. Original inputs are only ever read here, never modified."""
    vds1 = hl.vds.read_vds(hail_file_uri(args.input1))
    vds2 = hl.vds.read_vds(hail_file_uri(args.input2))

    ref_fields1 = vds_entry_field_order(vds1.reference_data)
    ref_fields2 = vds_entry_field_order(vds2.reference_data)
    var_fields1 = vds_entry_field_order(vds1.variant_data)
    var_fields2 = vds_entry_field_order(vds2.variant_data)

    if set(ref_fields1) != set(ref_fields2):
        raise ValueError(
            "input1 and input2 reference_data have different entry fields (not just different "
            f"order): {sorted(ref_fields1)} vs {sorted(ref_fields2)}"
        )
    if set(var_fields1) != set(var_fields2):
        raise ValueError(
            "input1 and input2 variant_data have different entry fields (not just different "
            f"order): {sorted(var_fields1)} vs {sorted(var_fields2)}"
        )

    ref_order = sorted(ref_fields1)
    var_order = sorted(var_fields1)

    ref_types1 = entry_field_types(vds1.reference_data)
    ref_types2 = entry_field_types(vds2.reference_data)
    var_types1 = entry_field_types(vds1.variant_data)
    var_types2 = entry_field_types(vds2.variant_data)

    canonical_ref_types = {f: resolve_canonical_type(f, ref_types1[f], ref_types2[f]) for f in ref_order}
    canonical_var_types = {f: resolve_canonical_type(f, var_types1[f], var_types2[f]) for f in var_order}

    aligned_dir = args.temp_dir / "aligned"
    input1_path = hail_file_uri(args.input1)
    input2_path = hail_file_uri(args.input2)

    needs_align1 = (
        ref_fields1 != ref_order
        or var_fields1 != var_order
        or any(ref_types1[f] != canonical_ref_types[f] for f in ref_order)
        or any(var_types1[f] != canonical_var_types[f] for f in var_order)
    )
    needs_align2 = (
        ref_fields2 != ref_order
        or var_fields2 != var_order
        or any(ref_types2[f] != canonical_ref_types[f] for f in ref_order)
        or any(var_types2[f] != canonical_var_types[f] for f in var_order)
    )

    if needs_align1:
        aligned1 = aligned_dir / "input1.vds"
        LOGGER.info("=== input1 schema differs from canonical (order and/or types); rewriting to %s", aligned1)
        remove_path(aligned1)
        aligned_vds1 = hl.vds.VariantDataset(
            reference_data=align_matrix_table_schema(vds1.reference_data, ref_order, canonical_ref_types),
            variant_data=align_matrix_table_schema(vds1.variant_data, var_order, canonical_var_types),
        )
        aligned_vds1.write(hail_file_uri(aligned1), overwrite=True)
        input1_path = hail_file_uri(aligned1)

    if needs_align2:
        aligned2 = aligned_dir / "input2.vds"
        LOGGER.info("=== input2 schema differs from canonical (order and/or types); rewriting to %s", aligned2)
        remove_path(aligned2)
        aligned_vds2 = hl.vds.VariantDataset(
            reference_data=align_matrix_table_schema(vds2.reference_data, ref_order, canonical_ref_types),
            variant_data=align_matrix_table_schema(vds2.variant_data, var_order, canonical_var_types),
        )
        aligned_vds2.write(hail_file_uri(aligned2), overwrite=True)
        input2_path = hail_file_uri(aligned2)

    return input1_path, input2_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine two Hail VDS datasets into a single merged VDS.",
    )
    parser.add_argument("input1", type=Path, help="First input VDS path.")
    parser.add_argument("input2", type=Path, help="Second input VDS path.")
    parser.add_argument("output_vds", type=Path, help="Path to write the merged VDS.")
    parser.add_argument(
        "--temp-dir",
        type=Path,
        default=None,
        help="Temp/save directory for the combiner. Defaults to <output_vds>.tmp next to the output.",
    )
    parser.add_argument("--reference", default="GRCh38", help="Hail reference genome. Default: GRCh38.")
    parser.add_argument(
        "--call-fields",
        type=comma_separated_fields,
        default="",
        help=(
            "Colon-separated FORMAT call fields, passed to Hail's combiner call_fields option. "
            "Must match what the input VDS were built with. Default: none (matches Hail combiner default "
            "of no call fields for VDS-to-VDS merges built the same way)."
        ),
    )
    parser.add_argument(
        "--whole-genome",
        action="store_true",
        help="Use genome default intervals instead of exome default intervals.",
    )
    parser.add_argument(
        "-M",
        "--spark-memory",
        default="",
        help="Spark driver/executor memory (e.g. '8g'), used in PYSPARK_SUBMIT_ARGS. Empty to skip.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Remove an existing output_vds before writing.")
    parser.add_argument(
        "--verbose", action="store_true", help="Enable debug logging and print VDS stats when finished."
    )
    args = parser.parse_args(argv)

    if args.temp_dir is None:
        args.temp_dir = args.output_vds.with_suffix(args.output_vds.suffix + ".tmp")

    return args


def init_hail(args: argparse.Namespace) -> None:
    if args.spark_memory:
        os.environ["PYSPARK_SUBMIT_ARGS"] = (
            f"--driver-memory {args.spark_memory} --executor-memory {args.spark_memory} pyspark-shell"
        )
        LOGGER.info("=== PYSPARK_SUBMIT_ARGS=%s", os.environ["PYSPARK_SUBMIT_ARGS"])

    LOGGER.info("=== Initializing Hail (reference=%s, tmp_dir=%s)", args.reference, args.temp_dir)
    hl.init(tmp_dir=hail_file_uri(args.temp_dir))
    hl.default_reference(args.reference)


def validate_inputs(args: argparse.Namespace) -> None:
    for path, label in ((args.input1, "input1"), (args.input2, "input2")):
        if not path.exists():
            raise FileNotFoundError(f"{label} VDS does not exist: {path}")
    if args.input1.resolve() == args.input2.resolve():
        raise ValueError("input1 and input2 must be different VDS paths")
    if args.output_vds.resolve() in (args.input1.resolve(), args.input2.resolve()):
        raise ValueError(
            "output_vds must not be the same path as input1 or input2 "
            "(this would delete an input before it's read)"
        )

    validate_input_manifest(args.input1, args, "input1")
    validate_input_manifest(args.input2, args, "input2")

    if args.output_vds.exists():
        if not args.overwrite:
            raise FileExistsError(
                f"Output VDS already exists: {args.output_vds} (use --overwrite to replace it)"
            )
        LOGGER.info("=== Removing existing output because --overwrite was set: %s", args.output_vds)
        remove_path(args.output_vds)
        remove_path(vds_success_marker_path(args.output_vds))


def combine_vds(args: argparse.Namespace) -> None:
    validate_inputs(args)

    save_path = args.temp_dir / "combiner-state"
    remove_path(save_path)

    input1_path, input2_path = align_vds_for_combiner(args)
    vds_paths = [input1_path, input2_path]
    LOGGER.info("=== Combining VDS inputs %s into %s", vds_paths, args.output_vds)

    start = time.monotonic()
    combiner = hl.vds.new_combiner(
        output_path=hail_file_uri(args.output_vds),
        vds_paths=vds_paths,
        save_path=hail_file_uri(save_path),
        temp_path=hail_file_uri(args.temp_dir),
        use_genome_default_intervals=args.whole_genome,
        use_exome_default_intervals=(not args.whole_genome),
        call_fields=args.call_fields,
    )
    combiner.run()
    mark_vds_success(args.output_vds, final_manifest(args))
    LOGGER.info("=== Finished combining in %.1f seconds: %s", time.monotonic() - start, args.output_vds)
    LOGGER.info("=== Wrote success marker: %s", vds_success_marker_path(args.output_vds))

    if args.verbose:
        vds = hl.vds.read_vds(hail_file_uri(args.output_vds))
        print(vds.variant_data.count())  # (variants, samples)
        vds.variant_data.rows().describe()


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    configure_logging(args.verbose)
    LOGGER.info("=== combine_vds.py version: %s", SCRIPT_VERSION)
    init_hail(args)
    combine_vds(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
