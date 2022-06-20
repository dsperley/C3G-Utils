#!/bin/bash

## add ambiguous reads to uniquely mapping reads

for i in $(ls alignments_with_ambiguous_reads/*mr | cut -f2 -d "/" | cut -f1 -d ".")
do

	mv alignments_with_ambiguous_reads/${i}_amb_reads.mr alignments_with_ambiguous_reads/${i}.mr
	cat alignments_with_ambiguous_reads/${i}.mr alignments_with_ambiguous_reads/${i}_amb_reads.mr_1_ambiguous alignments_with_ambiguous_reads/${i}_amb_reads.mr_1_ambiguous > alignments_with_ambiguous_reads/${i}_amb_reads.mr

done	

	
