#!/bin/bash


sbatch --job-name=walt --account=rrg-bourqueg-ad --time=8:00:00 --mem-per-cpu=4000M --cpus-per-task=8 --mail-type=ALL --mail-user=danielle.perley@mcgill.ca  --output=walt_out_%A_%a.txt --array=1-12 ./Walt_alignment_amb_reads.sh | grep [0-9] | cut -f4 -d " "


## combined uniquely mapping and ambiguous reads
bash combine_reads.sh

j2_ids=$(sbatch --job-name=deduplicated --account=rrg-bourqueg-ad --time=10:00:00 --mem=250G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca --output=dedup_out_%A_%a.txt --array=1-12 ./deduplicate.sh alignments_with_ambiguous_reads | grep [0-9] | cut -f4 -d " " )


j3_ids=$(sbatch --job-name=methCounts  --account=rrg-bourqueg-ad  --time=8:00:00 --mem=16G  --mail-type=ALL --mail-user=danielle.perley@mcgill.ca --output=methCounts_out_%A_%a.txt --array=1-12 --dependency=afterok:$j2_ids ./methcounts.sh alignments_with_ambiguous_reads 
| grep [0-9] | cut -f4 -d " ")

## methylKit

module load mugqic/R_Bioconductor/4.1.0_3.13

j4_ids=$(sbatch --job-name=meth_load --account=rrg-bourqueg-ad --time=3:00:00 --mem=100G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --dependency=afterok:$j3_ids --output=meth_load_%j.txt ./Methylation_analysis_walt_create_methyl_obj.R Methylation_analysis_no_amb_reads "amb_reads_methylKit.txt" | grep [0-9] | cut -f4)

sbatch --job-name=meth_analysis --cpus-per-task=8 --account=rrg-bourqueg-ad --time=24:00:00 --dependency=afterok:$j4_ids --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_%j.txt ./Methylation_analysis_walt.R Methylation_analysis_no_amb_reads


#### for timepoint3 removed, just PCA analysis
sbatch --job-name=meth_analysis_amb_tp3_rm --cpus-per-task=8 --account=rrg-bourqueg-ad --time=4:00:00 --mem-per-cpu=14G --mail-type=ALL --mail-user=danielle.perley@mcgill.ca -n 1 --output=meth_out_amb_reads_tp3_rm_%j.txt ./Methylation_analysis_walt_3TP_rm.R  Methylation_analysis_amb_reads  Methylation_analysis_amb_reads_timepoint3_rm | grep [0-9] | cut -f4 -d " "





