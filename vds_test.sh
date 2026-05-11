#!/bin/bash
set -e
cpu=12
mem="12G"

export PYSPARK_SUBMIT_ARGS="--driver-memory $mem --executor-memory $mem pyspark-shell"
#bsubrun -n $cpu -M $mem
./vds-from-gvcf.py gVCFs_small/ vds_temp vds_final --whole-genome