#!/bin/bash

mkdir -p slurm_out/locatit ;
for i in `ls genpipes/alignment/*/*.sorted.realigned.bam`
do
	ID=`echo $i | cut -d/ -f3`

	echo "echo -e \"#!/bin/bash\nmodule purge && module load java/11.0.2 && \
java -Xmx24G -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Dsamjdk.buffer_size=4194304 -jar \
/lustre06/project/6061810/dsperley/bin/agent3.0/lib/locatit-2.0.5.jar \
-PM:xm,Q:xq,q:nQ,r:nR \
-q 20 \
-m 1 \
-c 2500 \
-S \
-IB \
-OB \
-C \
-i \
-L \
-b SureSelectHumanAllExonV7.Target.b38.bed \
-o genpipes/alignment/${ID}/${ID} \
$i raw_reads/${ID}/*MBC_0.txt.gz\" |  \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A rrg-vmooser --chdir=\`pwd\` -o slurm_out/locatit/${ID}.locit.out --job-name=${ID}.locit --mem=25G --time=24:00:00 -N 1 -n 6" 

done

