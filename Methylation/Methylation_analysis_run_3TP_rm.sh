#!/bin/bash

#SBATCH --job-name=meth_3TP_rm
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --output=meth_out_3TP_rm_%j.txt


module load mugqic/R_Bioconductor/4.1.0_3.13

Rscript Methylation_analysis_072821_TP3_rm.R 
