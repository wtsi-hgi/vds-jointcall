#!/usr/bin/env python3
"""Report missing values for every row and entry field of a Hail MatrixTable.

Intended for MatrixTables produced by ``vds-to-mt.py``. The field list is taken
from the MatrixTable schema, so it reflects whatever annotations the source
gVCFs and the export actually produced.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import time
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import hail as hl

LOGGER = logging.getLogger("validate-mt")

# Subsets of present entries reported for every entry field, in output order.
ENTRY_SUBSETS = ("all", "hom_ref", "non_ref", "gt_missing")
GT_SUBSETS = ENTRY_SUBSETS[1:]


@dataclass(slots=True, frozen=True)
class Config:
    input_mt: str
    reference: str
    tmp_dir: str
    spark_memory: str
    output_tsv: str | None


@dataclass(slots=True, frozen=True)
class Field:
    path: str
    dtype: str
    expr: hl.Expression
    is_key: bool = False


@dataclass(slots=True, frozen=True)
class MissingStat:
    level: str
    field: str
    dtype: str
    subset: str
    n_total: int
    n_missing: int

    @property
    def pct_missing(self) -> float | None:
        """Percentage of missing values, or None for an empty subset."""
        return 100.0 * self.n_missing / self.n_total if self.n_total else None


def format_pct(pct: float | None, digits: int = 2) -> str:
    return "NA" if pct is None else f"{pct:.{digits}f}"


def has_uri_scheme(path: str) -> bool:
    return urlparse(path).scheme != ""


def hail_path(path: str) -> str:
    if has_uri_scheme(path):
        return path
    return Path(path).resolve().as_uri()


def init_hail(config: Config) -> None:
    if config.spark_memory != "":
        os.environ["PYSPARK_SUBMIT_ARGS"] = (
            f"--driver-memory {config.spark_memory} --executor-memory {config.spark_memory} pyspark-shell"
        )
        LOGGER.info("=== PYSPARK_SUBMIT_ARGS=%s", os.environ["PYSPARK_SUBMIT_ARGS"])

    LOGGER.info("=== Initializing Hail with reference=%s tmp_dir=%s", config.reference, config.tmp_dir)
    hl.init(tmp_dir=hail_path(config.tmp_dir))
    hl.default_reference(config.reference)


def iter_leaf_fields(expr: hl.StructExpression, prefix: str = "", keys: Sequence[str] = ()) -> Iterator[Field]:
    """Yield non-struct fields, descending into nested structs (e.g. ``info.AC``)."""
    for name, sub in expr.items():
        path = f"{prefix}{name}"
        if isinstance(sub, hl.StructExpression) and len(sub) > 0:
            yield from iter_leaf_fields(sub, prefix=f"{path}.")
        else:
            yield Field(
                path=path, dtype=str(sub.dtype).replace(" ", ""), expr=sub, is_key=name in keys and prefix == ""
            )


def is_missing_expr(expr: hl.Expression) -> hl.BooleanExpression:
    """Missing value, or a container (array/set) with at least one missing element."""
    if isinstance(expr.dtype, hl.tarray | hl.tset):
        return hl.is_missing(expr) | hl.any(lambda x: hl.is_missing(x), expr)
    return hl.is_missing(expr)


def print_schema(row_fields: Sequence[Field], entry_fields: Sequence[Field]) -> None:
    print("\n=== Variant-level (row) fields")
    for f in row_fields:
        print(f"  {f.path:<30} {f.dtype}{'  [key]' if f.is_key else ''}")
    print("\n=== Genotype-level (entry) fields")
    for f in entry_fields:
        print(f"  {f.path:<30} {f.dtype}")


def row_missing_stats(mt: hl.MatrixTable, fields: Sequence[Field]) -> list[MissingStat]:
    # Field paths may contain dots, so aggregate under positional names.
    aggs = {f"f{i}": hl.agg.count_where(is_missing_expr(f.expr)) for i, f in enumerate(fields)}
    result = mt.aggregate_rows(hl.struct(n_rows=hl.agg.count(), **aggs))
    return [
        MissingStat("row", f.path, f.dtype, "all_rows", result.n_rows, result[f"f{i}"]) for i, f in enumerate(fields)
    ]


def entry_missing_stats(
    mt: hl.MatrixTable, fields: Sequence[Field], n_rows: int, n_cols: int
) -> tuple[MissingStat, list[MissingStat], bool]:
    """Return the missing-entry summary and per-field stats computed over present entries only."""
    # Hail stores missing entries (sparse cells, unfiltered entries) as a defined struct whose fields
    # are all missing, so hl.is_defined(mt.entry) cannot detect them. Treat an entry as present when
    # at least one of its top-level fields is defined.
    present = hl.any([hl.is_defined(mt.entry[name]) for name in mt.entry]) if len(mt.entry) > 0 else hl.bool(False)
    has_gt = "GT" in mt.entry and mt.GT.dtype == hl.tcall

    # GT classes are evaluated only on present entries; a present entry with missing GT is "gt_missing".
    subset_conditions: dict[str, hl.BooleanExpression] = {"all": present}
    if has_gt:
        subset_conditions["hom_ref"] = present & hl.coalesce(mt.GT.is_hom_ref(), False)
        subset_conditions["non_ref"] = present & hl.coalesce(mt.GT.is_non_ref(), False)
        subset_conditions["gt_missing"] = present & hl.is_missing(mt.GT)

    aggs: dict[str, hl.Expression] = {f"n_{name}": hl.agg.count_where(cond) for name, cond in subset_conditions.items()}
    for i, f in enumerate(fields):
        missing = is_missing_expr(f.expr)
        for name, cond in subset_conditions.items():
            aggs[f"f{i}_{name}"] = hl.agg.count_where(cond & missing)

    result = mt.aggregate_entries(hl.struct(**aggs))

    # aggregate_entries skips filtered entries, so take the total from the dimensions;
    # filtered entries count as missing entries.
    n_cells = n_rows * n_cols
    summary = MissingStat("entry", "*", "*", "entries", n_cells, n_cells - result.n_all)
    totals = {name: result[f"n_{name}"] for name in subset_conditions}

    stats = [
        MissingStat("entry", f.path, f.dtype, name, totals[name], result[f"f{i}_{name}"])
        for i, f in enumerate(fields)
        for name in subset_conditions
    ]
    return summary, stats, has_gt


def print_table(header: Sequence[str], rows: Sequence[Sequence[str]], numeric_from: int) -> None:
    """Print columns separated by tabs, padded with spaces; columns from ``numeric_from`` are right-aligned."""
    widths = [max(len(row[i]) for row in [header, *rows]) for i in range(len(header))]

    def fmt(row: Sequence[str]) -> str:
        cells = [
            cell.rjust(width) if i >= numeric_from else cell.ljust(width)
            for i, (cell, width) in enumerate(zip(row, widths, strict=True))
        ]
        return "\t".join(cells).rstrip()

    print(fmt(header))
    for row in rows:
        print(fmt(row))


def print_entry_summary(summary: MissingStat, stats: Sequence[MissingStat]) -> None:
    totals = {s.subset: s.n_total for s in stats}
    n_present = totals.get("all", 0)
    print("\n=== Entries")
    print(f"Entries: total={summary.n_total} missing={summary.n_missing} missing_pct={format_pct(summary.pct_missing)}")
    parts = [f"present={n_present}"]
    for subset in GT_SUBSETS:
        if subset in totals:
            pct = 100.0 * totals[subset] / n_present if n_present else None
            parts.append(f"{subset}={totals[subset]} {subset}_pct={format_pct(pct)}")
    print("Present entries: " + " ".join(parts))


def print_row_report(stats: Sequence[MissingStat], n_rows: int) -> None:
    print(f"\n=== Missing values in variant-level fields (n_variants={n_rows})")
    rows = [[s.field, s.dtype, str(s.n_missing), format_pct(s.pct_missing)] for s in stats]
    print_table(["field", "type", "missing_n", "missing_pct"], rows, numeric_from=2)


def print_entry_report(stats: Sequence[MissingStat], has_gt: bool) -> None:
    subsets = ENTRY_SUBSETS if has_gt else ENTRY_SUBSETS[:1]
    by_field: dict[str, dict[str, MissingStat]] = {}
    for s in stats:
        by_field.setdefault(s.field, {})[s.subset] = s

    print("\n=== Missing values in genotype-level fields (present entries only)")
    print(
        "# all: all present entries; hom_ref / non_ref / gt_missing: present entries by GT; pct is NA for empty subsets"
    )
    if not has_gt:
        print("# no GT call field found; hom-ref / non-ref breakdown skipped")
    header = ["field", "type"] + [f"{name}_{col}" for name in subsets for col in ("n", "pct")]
    rows = []
    for field, per_subset in by_field.items():
        row = [field, per_subset[subsets[0]].dtype]
        for name in subsets:
            stat = per_subset[name]
            row += [str(stat.n_missing), format_pct(stat.pct_missing)]
        rows.append(row)
    print_table(header, rows, numeric_from=2)


def write_tsv(path: str, stats: Sequence[MissingStat]) -> None:
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["level", "field", "type", "subset", "n_total", "n_missing", "pct_missing"])
        for s in stats:
            writer.writerow([s.level, s.field, s.dtype, s.subset, s.n_total, s.n_missing, format_pct(s.pct_missing, 4)])
    LOGGER.info("=== Wrote TSV report: %s", path)


def validate_matrix_table(config: Config) -> None:
    start = time.monotonic()
    LOGGER.info("=== Reading MatrixTable: %s", config.input_mt)
    mt = hl.read_matrix_table(hail_path(config.input_mt))

    n_rows, n_cols = mt.count()
    print(f"\nMatrixTable: {config.input_mt}")
    print(f"  variants (rows): {n_rows:,}")
    print(f"  samples (cols):  {n_cols:,}")

    row_fields = list(iter_leaf_fields(mt.row, keys=list(mt.row_key)))
    entry_fields = list(iter_leaf_fields(mt.entry))
    print_schema(row_fields, entry_fields)

    LOGGER.info("=== Counting missing values in %d row field(s)", len(row_fields))
    row_stats = row_missing_stats(mt, row_fields)

    LOGGER.info("=== Counting missing values in %d entry field(s)", len(entry_fields))
    entry_summary, entry_stats, has_gt = entry_missing_stats(mt, entry_fields, n_rows, n_cols)

    print_row_report(row_stats, n_rows)
    print_entry_summary(entry_summary, entry_stats)
    print_entry_report(entry_stats, has_gt)

    if config.output_tsv is not None:
        write_tsv(config.output_tsv, [*row_stats, entry_summary, *entry_stats])
    LOGGER.info("Finished validating %s in %.1f seconds", config.input_mt, time.monotonic() - start)


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_args(argv: Sequence[str] | None = None) -> Config:
    parser = argparse.ArgumentParser(
        description=(
            "Report missing values for every variant-level and genotype-level field of a Hail MatrixTable "
            "(e.g. one produced by vds-to-mt.py). Array/set fields such as PL count as missing if any element "
            "is missing."
        ),
    )
    parser.add_argument("input_mt", help="Input Hail MatrixTable path or URI.")
    parser.add_argument("--reference", default="GRCh38", help="Hail reference genome. Default: GRCh38.")
    parser.add_argument(
        "--tmp-dir", default="hail-tmp", help="Hail/Spark temporary directory path or URI. Default: hail-tmp."
    )
    parser.add_argument(
        "-M",
        "--spark-memory",
        type=str,
        default="",
        help="Spark driver and executor memory used in PYSPARK_SUBMIT_ARGS. Use an empty value to skip setting it.",
    )
    parser.add_argument(
        "--output-tsv",
        default=None,
        help="Optional local path for a long-format TSV with all missing-value statistics.",
    )
    args = parser.parse_args(argv)
    return Config(
        input_mt=args.input_mt,
        reference=args.reference,
        tmp_dir=args.tmp_dir,
        spark_memory=args.spark_memory,
        output_tsv=args.output_tsv,
    )


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    config = parse_args(argv)
    init_hail(config)
    validate_matrix_table(config)
    return 0


if __name__ == "__main__":
    sys.exit(main())
