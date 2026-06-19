#!/bin/bash
set -e
mem="12G"

python vds-from-gvcf.py gVCFs_small vds_small --whole-genome --shard-size 4
