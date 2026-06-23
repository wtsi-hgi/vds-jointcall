# VDS-Jointcall

VDS-jointcall is a tool that makes joint calling of WGS/WES samples using Hail VDS format.

It is designed to work on large-scale datasets (up to 10_000 full genomes),
designed to accept gVCFs called by DeepVariant.

It required the Hail library.

## Export a VDS to MatrixTable

Use `vds-to-mt.py` to convert a Hail VDS into a filtered sparse Hail MatrixTable:

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt --overwrite
```

By default, the exporter:

- keeps variants with 10 total alleles or fewer
- adds `AD` and `GT` entry fields from VDS local allele fields
- creates a merged sparse MatrixTable
- splits multi-allelic sites
- keeps sites with at least one non-reference genotype
- drops the `gvcf_info` row field

The reference sequence FASTA is configurable when `split_multi_hts` needs it:

```bash
python3 vds-to-mt.py path/to/input.vds path/to/output.mt \
  --reference-sequence file:///path/to/hs38DH.fa \
  --max-alleles 10 \
  --overwrite
```
