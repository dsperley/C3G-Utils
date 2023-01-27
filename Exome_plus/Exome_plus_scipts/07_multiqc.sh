#!/bin/bash

#SBATCH -A rrg-vmooser
#SBATCH --time 1:00:00
#SBATCH --mem 2G



#export RAP_ID=vmooser
module load mugqic/python/3.7.3 mugqic/MultiQC/1.9

multiqc -o genpipes/metrics  genpipes/metrics/dna
