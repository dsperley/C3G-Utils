#!/bin/bash

#SBATCH --job-name=mapper
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danielle.perley@mcgill.ca
#SBATCH -o /home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401/logs/%j



### better to run in mapper directory
module load mugqic/mirdeep2/0.0.8


PROJECT="/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401"
TRIMMED=${PROJECT}/02_Trim
MAPPED=${PROJECT}/04_Mapped
LOGS=${PROJECT}/logs

if [ ! -d $MAPPED ]; then mkdir $MAPPED; fi
if [ ! -d $LOGS ]; then mkdir $LOGS; fi


#for i in $TRIMMED/*fastq.gz
#do
#
#	name=$(basename $i .gz)
#	gzip -c -d $i > $HOME/scratch/$name
#
#	## add location to config file
#	sample=$(echo $name | cut -f2 -d_| cut -f1 -d.)
#	echo -e "$HOME/scratch/$name\t${sample}" >> $PROJECT/config.txt
#	
#done


### -e fastq input                                                                                                                     ### -o threads for bowtie                                                                                                             ### -h parse to fasta format                     
### -m collapse reads                                                                                                                 
### -p genome index                                                                                                                    

mapper.pl $PROJECT/config.txt -d -e -o 6 -h -m -p ${PROJECT}/Ref/CamFam3.1_Canid_virus -s ${MAPPED}/reads.fa -t ${MAPPED}/reads_vs_genome.arf -v 2> $LOGS/mapper.out
