#!/bin/bash
set -e
mem="12G"

./vds-from-gvcf.py gVCFs_small/ vds_final_small --whole-genome --spark-memory "$mem"
