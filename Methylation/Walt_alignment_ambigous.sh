#!/bin/bash

#SBATCH --job-name=walt_amb
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000M
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH --output=walt_out_%A_%a.txt
#SBATCH --array=1-12


SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

## should have saved the fastq files in $HOME/scratch
#gunzip -cd trim/${SAMP}.trim.pair1.fastq.gz > trim/${SAMP}.trim.pair1.fastq
#gunzip -cd trim/${SAMP}.trim.pair2.fastq.gz > trim/${SAMP}.trim.pair2.fastq


/home/dsperley/walt/bin/walt -a -i ../../genomes/Okis_V2_genomic.dbindex -1 trim/${SAMP}.trim.pair1.fastq -2 trim/${SAMP}.trim.pair2.fastq -m 6 -k 3 -N 5000000 -o alignments_with_ambiguous_reads/${SAMP}_amb.mr -t 8
