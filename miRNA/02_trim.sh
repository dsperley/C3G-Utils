#!/bin/bash

#SBATCH --job-name=cutadapt
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=20:00
#SBATCH --mem-per-cpu=400M
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH --array=1-15
#SBATCH -o ./logs/%A_%a.out

module unload mugqic/python/2.7.14   
module load python/3.9.6 mugqic/cutadapt/2.10 fastqc/0.11.9


PROJECT="/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401"
SHEET="${PROJECT}/Adapters.txt"
OUT="${PROJECT}/02_Trim"

if [ ! -d "$OUT" ]; then mkdir "$OUT"; fi


adapter=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SHEET} | cut -f 2)
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SHEET} | cut -f 1)

cutadapt -j 5 -a $adapter -m 18 -M 30 -q 20 -o $OUT/${sample}.trim.fastq.gz ${PROJECT}/fastq/${sample}_R1.fastq.gz > $OUT/${sample}.trim.report.txt

fastqc -o ${OUT} $OUT/${sample}.trim.fastq.gz
