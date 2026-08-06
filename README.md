# VDS-Jointcall

VDS-jointcall is a set of scripts for joint calling of WGS/WES samples using the
[Hail VDS](https://hail.is/docs/0.2/vds/index.html) format.
It is designed to work on large-scale datasets (up to 10 000 full genomes) and
accepts gVCFs produced by DeepVariant.

The typical workflow is:

1. **Combine gVCFs into a VDS** — `vds-from-gvcf.py`
2. **Convert the VDS to a MatrixTable** — `vds-to-mt.py`
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

Use `vds-to-mt.py` to convert a Hail VDS into a filtered sparse Hail MatrixTable:

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt --overwrite
```

By default, the exporter:

- Keeps variants with 10 total alleles or fewer
- Adds `AD` and `GT` entry fields from VDS local allele fields
- Creates a merged sparse MatrixTable
- Splits multi-allelic sites
- Keeps sites with at least one non-reference genotype
- Drops the `gvcf_info` row field

The reference sequence FASTA is configurable when `split_multi_hts` needs it:

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt \
  --reference-sequence file:///path/to/hs38DH.fa \
  --max-alleles 10 \
  --overwrite
```

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
