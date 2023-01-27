#!/bin/bash

export RAP_ID=rrg-vmooser

mkdir -p slurm_out/agent_trim

for i in `ls -d raw_reads/*/*R1_001.fastq.gz`
do
	R2=`echo $i | sed -e 's#R1#R2#g'`
	ID=`echo $i | cut -d/ -f1-2`
	SID=`echo $i | cut -d/ -f2 | sed 's/Sample_//`

	echo "sbatch --mail-user=$JOB_MAIL --mail-type=end,fail --account=$RAP_ID --job-name=\"agent_trim.${SID}\" -o slurm_out/agent_trim/${SID}.trim.%j.out -n 1 -c 5 --mem=25G --time=6:00:00 --wrap=\"\
	module purge && \
	module load java/11.0.2 && \
	java -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Dsamjdk.buffer_size=4194304 -Xmx24G -jar \
	/lustre06/project/6061810/dsperley/bin/agent3.0/lib/trimmer-3.0.3.jar \
	-v2 -IDEE_FIXE \
	-fq1 $i \
	-fq2 $R2 \
	-out_loc $ID\""
done
