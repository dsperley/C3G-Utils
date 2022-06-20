#!/bin/bash

#SBATCH --job-name=mirdeep
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH -o /home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401/logs/%j

module load mugqic/mirdeep2/0.0.8
PROJECT=/project/6007512/C3G/projects/Pearson_canine_miRNA_R001401
OUT=$PROJECT/05_miRNA_detection_quantification
LOGS=${PROJECT}/logs

if [ ! -d $OUT ]; then mkdir $OUT; fi

miRDeep2.pl $PROJECT/04_Mapped/reads.fa $PROJECT/Ref/Canis_familiaris.CanFam3.1_Canid_herpesvirus_no_white_spaces.fa $PROJECT/04_Mapped/reads_vs_genome.arf $PROJECT/Ref/CanFam3.1_mature_no_white_space.fa none $PROJECT/Ref/CanFam3.1_hairpins_no_white_space.fa -r "$OUT/" 2 > $LOGS/mirdeep.out 
