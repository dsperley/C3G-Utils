#!/bin/bash

#SBATCH --job-name=walt
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000M
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH --output=walt_out_%A_%a.txt
#SBATCH --array=1-12


SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

## moved these files to $HOME/scratch
#gunzip -cd trim/${SAMP}.trim.pair1.fastq.gz > trim/${SAMP}.trim.pair1.fastq
#gunzip -cd trim/${SAMP}.trim.pair2.fastq.gz > trim/${SAMP}.trim.pair2.fastq


/home/dsperley/walt/bin/walt -i ../../genomes/Okis_V2_genomic.dbindex -1 $HOME/scratch/${SAMP}.trim.pair1.fastq -2 $HOME/scratch/${SAMP}.trim.pair2.fastq -m 6 -k 3 -N 5000000 -o alignments/${SAMP}.mr -t 8
