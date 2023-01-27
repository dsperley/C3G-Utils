#!/bin/bash

export RAP_ID=rrg-vmooser

mkdir -p slurm_out/cleanSam ;
for i in `ls genpipes/alignment/*/*.sorted.dup.recal.bam`
do 
 
	#ID=`echo $i | cut -d/ -f1-2` 
	SID=`echo $i | cut -d/ -f3` 
	
	out=${i%.*}.cleaned.bam
	echo "sbatch --mail-user=$JOB_MAIL --mail-type=end,fail -o slurm_out/cleanSam/${SID}_%j.out --account=$RAP_ID --job-name=\"clean.sam.${SID}\" -n 1 -c 10 --mem=36G --time=6:00:00 --wrap=\"\
	module purge && \
	module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/4.1.8.1 && \
	gatk --java-options \\\"-Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx36G\\\" \
	CleanSam \
	 -I $i \
	 -O $out \
	 -R /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
	 --CREATE_INDEX true
	\""
done
