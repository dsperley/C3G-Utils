#!/bin/bash

export RAP_ID=rrg-vmooser

mkdir -p slurm_out/cram

for i in `ls genpipes/alignment/*/*cleaned.bam`
do

	#ID=`echo $i | cut -d/ -f1-2`
	SID=`echo $i | cut -d/ -f3`

	echo "sbatch --mail-user=$JOB_MAIL --mail-type=end,fail --account=$RAP_ID --job-name=\"cram.${SID}\" -o slurm_out/cram/${SID}.cram.%j.out -n 1 -c 3 --mem=12G --time=3:00:00 --wrap=\"\
	module purge && \
	module load mugqic/samtools/1.14 && \
	samtools view -C --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa --threads 3 -o ${i%.*}.cram $i\""
done
