#!/bin/bash

path=/lustre06/project/6061810/dsperley/E100043100
for i in $(ls genpipes/alignment/*/*cleaned.bam)
do

	bai=$(echo $i | sed 's/bam/bai/')
	bam2=$(echo $i | sed 's/cleaned.//')

	 ln -sf ${path}/$i ${path}/${bam2}
	 ln -sf ${path}/$bai ${path}/${bam2}.bai

done

