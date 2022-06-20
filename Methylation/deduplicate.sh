#!/bin/bash
  
ALIGN_DIR="$1"

SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

BASE=$(basename $ALIGN_DIR/${SAMP}_amb_reads.mr .mr)

LC_ALL=C sort -k 1,1 -k2,2n -k3,3n -k6,6 -o $HOME/scratch/${BASE}_sorted.mr $ALIGN_DIR/${BASE}.mr

~/methpipe-4.1.1/build/duplicate-remover -S $ALIGN_DIR/${BASE}_dedup.stats.txt -o $ALIGN_DIR/${BASE}_dup_rm_sorted.mr $HOME/scratch/${BASE}_sorted.mr

