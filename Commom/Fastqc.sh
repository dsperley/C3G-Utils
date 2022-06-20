#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2375M
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH -n 1
#SBATCH --output=fastq_out_%j.txt



module unload mugqic/python/2.7.14
module load  python/3.9.6 mugqic/java/openjdk-jdk1.8.0_72 mugqic/fastqc/0.11.5 mugqic/MultiQC/1.10.1

OUT="01_FASTQC"

if [ ! -d "$OUT" ]; then mkdir "$OUT"; fi

fastqc -t 16 -o $OUT/raw_fastq/NS*R[12].fastq.gz

multiqc ${OUT}/*zip
