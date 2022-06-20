#!/bin/bash

PROJECT=/project/6007512/C3G/projects/Pearson_canine_miRNA_R001401
TRIM=$PROJECT/02_Trim
REPORT=$PROJECT/Reports

if [ ! -d $REPORT ]; then mkdir $REPORT; fi

### code to collect all the stats

echo -e "Sample\tRaw_Reads\tReads_Retained\tPercentage_Retained" > $REPORT/trim_summary.txt 
for i in $TRIM/*report.txt
do

	line1=$(grep "Total reads processed:" $i)
	total=$(echo $line1 | grep -E -o "[0-9]{2},[0-9]{3},[0-9]{3}"| sed 's/,//g')

	line2=$(grep "Reads written (passing filters)" $i)
	read_num=$(echo $line2 | grep -E -o "[0-9]{2},[0-9]{3},[0-9]{3}"| sed 's/,//g')
	perc=$(echo $line2 | grep -E -o "[0-9]{2}.[0-9]{1,3}%")
	sample=$(basename $i .trim.report.txt)

	echo -e "${sample}\t${total}\t${read_num}\t${perc}" >> $REPORT/trim_summary.txt

done	

#module load mugqic/R_Bioconductor/4.1.0_3.13

#Rscript $PROJECT/Scripts/alignment_stats.R $PROJECT/04_Mapped/reads_vs_genome.arf $PROJECT/Reports/alignment_stats.csv


### R code to combine
 trim <- read.delim("trim_summary.txt") 
 align <- read.delim("alignment.txt")   

 trim <- mutate(trim,sample=str_match(Sample, ".*(S\\d{2})\\.M.*")[,2])


