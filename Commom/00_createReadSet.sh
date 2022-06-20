## Using sheet

PROJECT="/project/6007512/C3G/projects/PKR_Cachexia_RNA-Seq_Vahab_R004281"
ADAPTER="CTGTCTCTTATACACATCT"
RUN_TYPE="SINGLE_END"


echo -e "Sample\tReadset\tLibrary\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM"> ../ReadSet.txt
while read line
do

	## single end
	sample=$(echo $i | cut -f1 -d "_")
	library=$(echo $i | cut -f1 -d "_")
	lane=$(echo $i | cut -f4 -d "_")
	run=$(echo $i | cut -f2 -d "_" | sed 's/A.*//')

	echo -e "${sample}\t${sample}\t${library}\t${RUN_TYPE}\t${run}\t${lane}\t${ADAPTER}\t\t33\t\t${PROJECT}/fastq/${i}\t\t" >> ../ReadSet.txt

done
	
