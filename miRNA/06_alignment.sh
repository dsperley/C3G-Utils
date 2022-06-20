#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH -o /home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401/logs/%A_%a
#SBATCH --array=1-7

## to generate bam files
module load mugqic/bowtie/1.2.2 mugqic/samtools/1.12 mugqic/python/3.7.3

PROJECT=/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401
#fq=$(sed -n ${SLURM_ARRAY_TASK_ID}p $PROJECT/config.txt | cut -f1)
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $PROJECT/files_to_sort.txt)
fq="$HOME/scratch/${sample}.trim.fastq"

#sample=$(basename $fq .trim.fastq)
#fq=$(sed -n ${SLURM_ARRAY_TASK_ID}p $PROJECT/files_to_sort.txt)
bowtie -p 6 -n 0 -e 80 -l 18 -a -m 5 -S --best --strata /home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401/Ref/CamFam3.1_Canid_virus $fq 2>$PROJECT/04_Mapped/${sample}_bowtie.log | samtools view -bS -F4 - > $PROJECT/04_Mapped/${sample}.bam 

#sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $PROJECT/files_to_sort.txt)
samtools sort -o $PROJECT/04_Mapped/${sample}.sort.bam $PROJECT/04_Mapped/${sample}.bam
samtools index $PROJECT/04_Mapped/${sample}.sort.bam 
