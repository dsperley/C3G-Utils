#!/bin/bash

module load mugqic/R_Bioconductor/4.1.0_3.13

## for analysis with no ambiguous reads

#step1_ids=$(sbatch --job-name=meth_load --account=rrg-bourqueg-ad --time=7:00:00 --mem=400G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_load_%j.txt ./Methylation_analysis_walt_create_methyl_obj.R Methylation_analysis_no_amb_reads "\\d_methylKit.txt" | grep [0-9] | cut -f4)

#sbatch --job-name=meth_analysis --cpus-per-task=8 --account=rrg-bourqueg-ad --time=16:00:00 --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_%j.txt ./Methylation_analysis_walt.R Methylation_analysis_no_amb_reads | grep [0-9] | cut -f4

sbatch --job-name=meth_analysis_no_amb_tp3_rm --cpus-per-task=8 --account=rrg-bourqueg-ad --time=4:00:00 --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_%j.txt ./Methylation_analysis_walt_3TP_rm.R Methylation_analysis_no_amb_reads  Methylation_analysis_no_amb_reads_timepoint3_rm | grep [0-9] | cut -f4 -d " "


## for analysis with ambiguous  reads

#step1a_ids=$(sbatch --job-name=meth_load_amb_reads --account=rrg-bourqueg-ad --time=7:00:00 --mem=400G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_load_amb_reads_%j.txt ./Methylation_analysis_walt_create_methyl_obj.R Methylation_analysis_amb_reads _amb_reads_methylKit.txt | grep [0-9] | cut -f4)

#sbatch --job-name=meth_analysis_amb_reads --cpus-per-task=8 --account=rrg-bourqueg-ad --time=16:00:00 --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_amb_reads_%j.txt ./Methylation_analysis_walt.R Methylation_analysis_amb_reads | grep [0-9] | cut -f4
  
sbatch --job-name=meth_analysis_amb_tp3_rm --cpus-per-task=8 --account=rrg-bourqueg-ad --time=4:00:00 --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_amb_reads_tp3_rm_%j.txt ./Methylation_analysis_walt_3TP_rm.R  Methylation_analysis_amb_reads  Methylation_analysis_amb_reads_timepoint3_rm | grep [0-9] | cut -f4 -d " "
  
