# VDS-Jointcall

VDS-jointcall is a set of scripts for joint calling of WGS/WES samples using the
[Hail VDS](https://hail.is/docs/0.2/vds/index.html) format.
It is designed to work on large-scale datasets (up to 10 000 full genomes) and
accepts gVCFs produced by DeepVariant.

The typical workflow is:

1. **Combine gVCFs into a VDS** — `vds-from-gvcf.py`
2. **Convert the VDS to an annotated MatrixTable** — `vds-to-mt.py`
3. (Optional) **Export a single sample from the MatrixTable to VCF** — `mt-to-vcf-1sample.py`

It requires the [Hail](https://hail.is/) library 
(see `pyproject.toml` for the exact version).

## Combine gVCFs into a VDS

Use `vds-from-gvcf.py` to combine all gVCF files in a local folder into a single
Hail VDS:

```bash
python3 vds-from-gvcf.py path/to/gvcf_dir path/to/output.vds
```

The script:

- discovers every `.g.vcf.gz` / `.g.vcf.bgz` file in the input folder
- reads the sample ID from each gVCF header
- groups gVCFs into shards (default 50 per shard) and combines each shard into a temporary VDS
- merges all shard VDS files into the final output VDS
- tracks progress with per-shard manifests so interrupted runs can be resumed

### Appending new samples to an existing VDS

To add new gVCFs to an already-combined VDS without rebuilding from scratch, use
`--append-to`:

```bash
python3 vds-from-gvcf.py path/to/new_gvcfs path/to/updated_output.vds \
  --append-to path/to/existing.vds
```

Only gVCFs whose sample IDs are not already present in the existing VDS will be
processed.

### Common options

| Option                     | Description                                          |
|----------------------------|------------------------------------------------------|
| `--shard-size N`           | Number of gVCFs per shard (default 50)               |
| `--whole-genome`           | Use whole-genome intervals instead of exome defaults |
| `--recursive`              | Search the gVCF folder recursively                   |
| `--overwrite`              | Delete existing outputs before writing               |
| `--temp-vds-dir DIR`       | Custom temporary folder for intermediate VDS files   |
| `-M`, `--spark-memory MEM` | Spark driver/executor memory (e.g. `16g`)            |

## Export a VDS to MatrixTable

Use `vds-to-mt.py` to convert a Hail VDS into a filtered, annotated dense Hail MatrixTable:

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt --overwrite
```

By default, the exporter:

- Keeps variants with 10 total alleles or fewer (`--max-alleles`)
- Keeps sites with at least one non-reference genotype
- Creates a dense MatrixTable (`hl.vds.to_dense_mt`): samples covered by a reference block at a
  variant site get a hom-ref (`0/0`) entry, so AN/AF and other per-site aggregates are correct
- Keeps multi-allelic sites as they are (no splitting)
- Builds entry fields (`GT`, `AD`, `PL`, ...) from the VDS local allele fields and site-level
  `info` fields, as `annotate_vds.py` does
- Drops the VDS local fields (`LA`, `LGT`, `LAD`, `LPL`, `LPGT`) and the `gvcf_info` entry field
  (keep it with `--keep-gvcf-info`)

The output is multi-allelic. Allele-specific `info` arrays (`AS_*`) and `info.AC`/`info.AF` have one
value per alternate allele (VCF `Number=A`). `bcftools norm -m-` splits them automatically; after
`hl.split_multi_hts` use `info.AS_X[a_index - 1]`, because Hail does not re-index row `info` arrays.
Densifying scans all reference blocks, so expect a longer runtime than a sparse export.

`--reference-sequence` registers a reference FASTA with Hail; the export itself does not need it:

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt \
  --reference-sequence file:///path/to/hs38DH.fa \
  --max-alleles 10 \
  --overwrite
```

### Annotations

`--annotations` selects what to build as a comma-separated list. Entry fields are plain names
(`GT`, `AD`), row `info` fields use the prefix `info.` (`info.QD`). `GT` is always built.

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt --annotations AD,DP,GQ,info.AC,info.AF,info.AN
```

- Default: all annotations below. Annotations whose source fields are missing in the VDS (common
  for DeepVariant gVCFs, which lack most GATK INFO fields) are skipped with a warning.
- An annotation named in `--annotations` fails the run if its source fields are missing.
  `--strict-annotations` makes the default list strict as well.
- `--no-annotations` builds only `GT`.
- `--list-annotations` prints all annotations with their required source fields and exits.
- `--fill-ref-dp` fills missing `DP` with `MIN_DP` (gnomAD convention); this changes `info.DP`.

Source fields are looked up as top-level VDS variant entry fields or inside `gvcf_info`.

| Annotation | Required sources | Notes |
|------------|------------------|-------|
| `GT` | `LGT`, `LA` | Always built |
| `AD` | `LAD`, `LA` | Global allele indices, missing for hom-ref entries |
| `PL` | `LPL`, `LA` | Global allele indices, missing for hom-ref entries |
| `PGT` | `LPGT`, `LA` | GATK only |
| `DP`, `GQ`, `MIN_DP`, `PID`, `PS`, `RGQ`, `SB` | same name | Copied from the entry, or from `gvcf_info` |
| `info.DP` | `DP` | Sum over all genotypes |
| `info.ReadPosRankSum`, `info.MQRankSum`, `info.BaseQRankSum` | same name | Median over non-ref genotypes |
| `info.MQ` | `RAW_MQandDP` | |
| `info.FS`, `info.SOR` | `SB` | Missing at sites where a non-ref `SB` is not in GATK notation (4 values) |
| `info.AS_MQ` | `AS_RAW_MQ`, `LAD`, `LA` | Per alt allele |
| `info.AS_ReadPosRankSum`, `info.AS_MQRankSum` | `AS_RAW_ReadPosRankSum` / `AS_RAW_MQRankSum`, `LA` | Per alt allele |
| `info.AS_QD` | `AS_QD`, `LA` | Per alt allele |
| `info.AS_FS`, `info.AS_SOR` | `AS_SB_TABLE`, `LA` | Per alt allele |
| `info.DB` | `DB` | |
| `info.QD` | `LPL`, `LA`, `LGT`, `DP` | Does not require the `PL` entry in the output |
| `info.AN`, `info.AC`, `info.AF` | `LGT`, `LA` | `AC`/`AF` per alt allele |

The site and allele-specific `info` fields are computed with the gnomAD helpers
`get_site_info_expr` and `get_as_info_expr`.

## Export a single sample to VCF

Use `mt-to-vcf-1sample.py` to extract one sample from a Hail MatrixTable and write it as a VCF file:

```bash
python3 mt-to-vcf-1sample.py path/to/input.mt path/to/output.vcf \
  --sample-id SAMPLE_NAME
```

The script filters the MatrixTable to the requested sample, drops the `gvcf_info` field, and exports the result to VCF.

This script is mostly used to extract a small VCF suitable for VEP annotation.
Therefore, it does not filter variations for non-ref genotypes and preserves all variants, 
even if in the exported sample they are in hom-ref.

### Common options

| Option                     | Description                                       |
|----------------------------|---------------------------------------------------|
| `--sample-id ID`           | Sample ID to extract (default: `NA18978`)         |
| `--num-partitions N`       | Number of partitions for the export (default 500) |
| `--overwrite`              | Overwrite an existing output VCF                  |
| `-M`, `--spark-memory MEM` | Spark driver/executor memory (e.g. `16g`)         |
