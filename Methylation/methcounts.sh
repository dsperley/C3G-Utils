#!/bin/bash


ALIGN_DIR="$1"
SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

BASE=$(basename $ALIGN_DIR/${SAMP}*dup_rm*.mr _dup_rm_sorted.mr)

if [ ! -d Methylation_calls ]; then mkdir Methylation_calls; fi

~/methpipe-4.1.1/build/methcounts -o Methylation_calls/${BASE}_ALL.meth.txt -c Okis_V2_genomic_names_fixed.fa  $ALIGN_DIR/${BASE}_dup_rm_sorted.mr


~/methpipe-4.1.1/build/symmetric-cpgs -o Methylation_calls/${BASE}_CpG.meth.txt Methylation_calls/${BASE}_ALL.meth.txt


awk -F $'\t' -v OFS=$'\t' '$6>0{$5=$5*100; $5=sprintf("%.2f",$5); print $1"."$2, $1,$2,"F",$6,$5,(100 - $5)}' Methylation_calls/${BASE}_CpG.meth.txt > Methylation_calls/${BASE}_methylKit.txt


