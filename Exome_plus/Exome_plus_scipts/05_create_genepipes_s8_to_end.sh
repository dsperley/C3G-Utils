#!/bin/bash

#export RAP_ID=vmooser

## when running pipeline, make sure to export RAP_ID before
module purge && module load mugqic/python/3.10.4 mugqic/genpipes/4.1.0 && dnaseq.py \
-c $MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.base.ini \
$MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.beluga.ini  \
$MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.exome.ini \
$MUGQIC_PIPELINES_HOME/pipelines/dnaseq/gatk4.ini \
$MUGQIC_PIPELINES_HOME/resources/genomes/config/Homo_sapiens.GRCh38.ini \
Exome_plus.ini -r readSet.txt -s 8,20-22 -o genpipes -g 05_dnaseq_s8_21_gatk4_d.sh

