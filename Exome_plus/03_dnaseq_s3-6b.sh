#!/bin/bash
# Exit immediately on error

set -eu -o pipefail

#-------------------------------------------------------------------------------
# DnaSeq SLURM Job Submission Bash script
# Version: 4.1.0
# Created on: 2022-08-04T11:31:47
# Steps:
#   bwa_mem_sambamba_sort_sam: 0 job... skipping
#   sambamba_merge_sam_files: 0 job... skipping
#   gatk_indel_realigner: 111 jobs
#   sambamba_merge_realigned: 0 job... skipping
#   TOTAL: 111 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/lustre06/project/6061810/dsperley/Run1/genpipes
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/DnaSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.1.0/pipelines/dnaseq/dnaseq.base.ini,/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.1.0/pipelines/dnaseq/dnaseq.exome.ini,/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.1.0/pipelines/dnaseq/dnaseq.beluga.ini,Exome_plus.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

#-------------------------------------------------------------------------------
# STEP: gatk_indel_realigner
#-------------------------------------------------------------------------------
STEP=gatk_indel_realigner
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_1_JOB_ID: gatk_indel_realigner.FD25671004_2-2130693
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671004_2-2130693
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671004_2-2130693.b1b5deb21ef5d57999d1277345d79c4f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671004_2-2130693.b1b5deb21ef5d57999d1277345d79c4f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671004_2-2130693/realign && \
touch alignment/FD25671004_2-2130693/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/FD25671004_2-2130693.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/FD25671004_2-2130693.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671004_2-2130693/FD25671004_2-2130693.sorted.realigned.bam
gatk_indel_realigner.FD25671004_2-2130693.b1b5deb21ef5d57999d1277345d79c4f.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_2_JOB_ID: gatk_indel_realigner.FD25671017_2-2134265
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671017_2-2134265
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671017_2-2134265.022526738664d961403d5c8fd07c294f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671017_2-2134265.022526738664d961403d5c8fd07c294f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671017_2-2134265/realign && \
touch alignment/FD25671017_2-2134265/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/FD25671017_2-2134265.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/FD25671017_2-2134265.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671017_2-2134265/FD25671017_2-2134265.sorted.realigned.bam
gatk_indel_realigner.FD25671017_2-2134265.022526738664d961403d5c8fd07c294f.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_3_JOB_ID: gatk_indel_realigner.FD25671018_2-2134267
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671018_2-2134267
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671018_2-2134267.4f5463093e680b21ab735d9442704c9a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671018_2-2134267.4f5463093e680b21ab735d9442704c9a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671018_2-2134267/realign && \
touch alignment/FD25671018_2-2134267/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/FD25671018_2-2134267.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/FD25671018_2-2134267.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671018_2-2134267/FD25671018_2-2134267.sorted.realigned.bam
gatk_indel_realigner.FD25671018_2-2134267.4f5463093e680b21ab735d9442704c9a.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_4_JOB_ID: gatk_indel_realigner.FD25671019_2-2134269
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671019_2-2134269
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671019_2-2134269.eb0b05fce087f4cac1cfc9f9f9d8898b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671019_2-2134269.eb0b05fce087f4cac1cfc9f9f9d8898b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671019_2-2134269/realign && \
touch alignment/FD25671019_2-2134269/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/FD25671019_2-2134269.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/FD25671019_2-2134269.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671019_2-2134269/FD25671019_2-2134269.sorted.realigned.bam
gatk_indel_realigner.FD25671019_2-2134269.eb0b05fce087f4cac1cfc9f9f9d8898b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_5_JOB_ID: gatk_indel_realigner.FD25671020_2-2134271
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671020_2-2134271
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671020_2-2134271.ec5babe3407054894944860a994c8afc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671020_2-2134271.ec5babe3407054894944860a994c8afc.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671020_2-2134271/realign && \
touch alignment/FD25671020_2-2134271/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/FD25671020_2-2134271.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/FD25671020_2-2134271.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671020_2-2134271/FD25671020_2-2134271.sorted.realigned.bam
gatk_indel_realigner.FD25671020_2-2134271.ec5babe3407054894944860a994c8afc.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_6_JOB_ID: gatk_indel_realigner.FD25671021_2-2134273
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671021_2-2134273
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671021_2-2134273.293e52f8cbaf22e5286c3a96493934da.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671021_2-2134273.293e52f8cbaf22e5286c3a96493934da.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671021_2-2134273/realign && \
touch alignment/FD25671021_2-2134273/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/FD25671021_2-2134273.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/FD25671021_2-2134273.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671021_2-2134273/FD25671021_2-2134273.sorted.realigned.bam
gatk_indel_realigner.FD25671021_2-2134273.293e52f8cbaf22e5286c3a96493934da.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_6_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_7_JOB_ID: gatk_indel_realigner.FD25671022_2-2134275
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671022_2-2134275
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671022_2-2134275.e3f71fef0bedb49dd5f0787d65ee91e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671022_2-2134275.e3f71fef0bedb49dd5f0787d65ee91e3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671022_2-2134275/realign && \
touch alignment/FD25671022_2-2134275/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/FD25671022_2-2134275.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/FD25671022_2-2134275.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671022_2-2134275/FD25671022_2-2134275.sorted.realigned.bam
gatk_indel_realigner.FD25671022_2-2134275.e3f71fef0bedb49dd5f0787d65ee91e3.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_7_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_7_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_8_JOB_ID: gatk_indel_realigner.FD25671023_2-2134277
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671023_2-2134277
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671023_2-2134277.68ffb8453cb7e49bc2c345b446f3576b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671023_2-2134277.68ffb8453cb7e49bc2c345b446f3576b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671023_2-2134277/realign && \
touch alignment/FD25671023_2-2134277/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/FD25671023_2-2134277.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/FD25671023_2-2134277.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671023_2-2134277/FD25671023_2-2134277.sorted.realigned.bam
gatk_indel_realigner.FD25671023_2-2134277.68ffb8453cb7e49bc2c345b446f3576b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_8_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_8_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_9_JOB_ID: gatk_indel_realigner.FD25671024_2-2134279
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671024_2-2134279
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671024_2-2134279.e6bb7d2a27f7b5746cd27a1ad856e4d1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671024_2-2134279.e6bb7d2a27f7b5746cd27a1ad856e4d1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671024_2-2134279/realign && \
touch alignment/FD25671024_2-2134279/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/FD25671024_2-2134279.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/FD25671024_2-2134279.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671024_2-2134279/FD25671024_2-2134279.sorted.realigned.bam
gatk_indel_realigner.FD25671024_2-2134279.e6bb7d2a27f7b5746cd27a1ad856e4d1.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_9_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_9_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_10_JOB_ID: gatk_indel_realigner.FD25671032_2-2130667
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25671032_2-2130667
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25671032_2-2130667.1ffac40e543c45f845349f413750fc97.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25671032_2-2130667.1ffac40e543c45f845349f413750fc97.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25671032_2-2130667/realign && \
touch alignment/FD25671032_2-2130667/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/FD25671032_2-2130667.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/FD25671032_2-2130667.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25671032_2-2130667/FD25671032_2-2130667.sorted.realigned.bam
gatk_indel_realigner.FD25671032_2-2130667.1ffac40e543c45f845349f413750fc97.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_10_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_10_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_11_JOB_ID: gatk_indel_realigner.FD25672250_2-2130755
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672250_2-2130755
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672250_2-2130755.d5416830dae0b7f12e27dc70f2e398ca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672250_2-2130755.d5416830dae0b7f12e27dc70f2e398ca.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672250_2-2130755/realign && \
touch alignment/FD25672250_2-2130755/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/FD25672250_2-2130755.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/FD25672250_2-2130755.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672250_2-2130755/FD25672250_2-2130755.sorted.realigned.bam
gatk_indel_realigner.FD25672250_2-2130755.d5416830dae0b7f12e27dc70f2e398ca.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_11_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_11_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_12_JOB_ID: gatk_indel_realigner.FD25672251_2-2130753
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672251_2-2130753
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672251_2-2130753.ad33fb20a25a97be6196db1c0867b1e1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672251_2-2130753.ad33fb20a25a97be6196db1c0867b1e1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672251_2-2130753/realign && \
touch alignment/FD25672251_2-2130753/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/FD25672251_2-2130753.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/FD25672251_2-2130753.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672251_2-2130753/FD25672251_2-2130753.sorted.realigned.bam
gatk_indel_realigner.FD25672251_2-2130753.ad33fb20a25a97be6196db1c0867b1e1.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_12_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_12_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_13_JOB_ID: gatk_indel_realigner.FD25672252_2-2130751
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672252_2-2130751
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672252_2-2130751.02677930979d91e66450e2560c06bb09.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672252_2-2130751.02677930979d91e66450e2560c06bb09.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672252_2-2130751/realign && \
touch alignment/FD25672252_2-2130751/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/FD25672252_2-2130751.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/FD25672252_2-2130751.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672252_2-2130751/FD25672252_2-2130751.sorted.realigned.bam
gatk_indel_realigner.FD25672252_2-2130751.02677930979d91e66450e2560c06bb09.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_13_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_13_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_14_JOB_ID: gatk_indel_realigner.FD25672253_2-2130699
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672253_2-2130699
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672253_2-2130699.23092e699eb8215c8f35df2de265b39f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672253_2-2130699.23092e699eb8215c8f35df2de265b39f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672253_2-2130699/realign && \
touch alignment/FD25672253_2-2130699/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/FD25672253_2-2130699.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/FD25672253_2-2130699.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672253_2-2130699/FD25672253_2-2130699.sorted.realigned.bam
gatk_indel_realigner.FD25672253_2-2130699.23092e699eb8215c8f35df2de265b39f.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_14_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_14_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_15_JOB_ID: gatk_indel_realigner.FD25672254_2-2130697
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672254_2-2130697
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672254_2-2130697.a2f8adb4a836b47f11151b5403fd066e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672254_2-2130697.a2f8adb4a836b47f11151b5403fd066e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672254_2-2130697/realign && \
touch alignment/FD25672254_2-2130697/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/FD25672254_2-2130697.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/FD25672254_2-2130697.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672254_2-2130697/FD25672254_2-2130697.sorted.realigned.bam
gatk_indel_realigner.FD25672254_2-2130697.a2f8adb4a836b47f11151b5403fd066e.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_15_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_15_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_16_JOB_ID: gatk_indel_realigner.FD25672256_2-2130694
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672256_2-2130694
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672256_2-2130694.7ad17238bb2557c4cc657a761db5a9f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672256_2-2130694.7ad17238bb2557c4cc657a761db5a9f0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672256_2-2130694/realign && \
touch alignment/FD25672256_2-2130694/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/FD25672256_2-2130694.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/FD25672256_2-2130694.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672256_2-2130694/FD25672256_2-2130694.sorted.realigned.bam
gatk_indel_realigner.FD25672256_2-2130694.7ad17238bb2557c4cc657a761db5a9f0.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_16_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_16_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_17_JOB_ID: gatk_indel_realigner.FD25672257_2-2130757
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672257_2-2130757
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672257_2-2130757.7ef7c6e8a7727732675b2675b74a4d8b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672257_2-2130757.7ef7c6e8a7727732675b2675b74a4d8b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672257_2-2130757/realign && \
touch alignment/FD25672257_2-2130757/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/FD25672257_2-2130757.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/FD25672257_2-2130757.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672257_2-2130757/FD25672257_2-2130757.sorted.realigned.bam
gatk_indel_realigner.FD25672257_2-2130757.7ef7c6e8a7727732675b2675b74a4d8b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_17_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_17_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_18_JOB_ID: gatk_indel_realigner.FD25672361_2-2130695
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672361_2-2130695
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672361_2-2130695.9296201fed34118b94adb145ce90481a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672361_2-2130695.9296201fed34118b94adb145ce90481a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672361_2-2130695/realign && \
touch alignment/FD25672361_2-2130695/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/FD25672361_2-2130695.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/FD25672361_2-2130695.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672361_2-2130695/FD25672361_2-2130695.sorted.realigned.bam
gatk_indel_realigner.FD25672361_2-2130695.9296201fed34118b94adb145ce90481a.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_18_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_18_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_19_JOB_ID: gatk_indel_realigner.FD25672362_2-2130696
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672362_2-2130696
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672362_2-2130696.6538476090afd4672d994c033293956b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672362_2-2130696.6538476090afd4672d994c033293956b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672362_2-2130696/realign && \
touch alignment/FD25672362_2-2130696/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/FD25672362_2-2130696.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/FD25672362_2-2130696.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672362_2-2130696/FD25672362_2-2130696.sorted.realigned.bam
gatk_indel_realigner.FD25672362_2-2130696.6538476090afd4672d994c033293956b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_19_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_19_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_20_JOB_ID: gatk_indel_realigner.FD25672363_2-2130698
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672363_2-2130698
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672363_2-2130698.d0c33fec892644969a87d9eac21e4ee9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672363_2-2130698.d0c33fec892644969a87d9eac21e4ee9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672363_2-2130698/realign && \
touch alignment/FD25672363_2-2130698/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/FD25672363_2-2130698.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/FD25672363_2-2130698.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672363_2-2130698/FD25672363_2-2130698.sorted.realigned.bam
gatk_indel_realigner.FD25672363_2-2130698.d0c33fec892644969a87d9eac21e4ee9.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_20_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_20_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_21_JOB_ID: gatk_indel_realigner.FD25672364_2-2130700
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672364_2-2130700
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672364_2-2130700.5c8168f39f62b94e08479abed8055507.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672364_2-2130700.5c8168f39f62b94e08479abed8055507.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672364_2-2130700/realign && \
touch alignment/FD25672364_2-2130700/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/FD25672364_2-2130700.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/FD25672364_2-2130700.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672364_2-2130700/FD25672364_2-2130700.sorted.realigned.bam
gatk_indel_realigner.FD25672364_2-2130700.5c8168f39f62b94e08479abed8055507.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_21_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_21_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_22_JOB_ID: gatk_indel_realigner.FD25672365_2-2130752
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672365_2-2130752
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672365_2-2130752.899f2e0a1d7988810045196bab741e70.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672365_2-2130752.899f2e0a1d7988810045196bab741e70.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672365_2-2130752/realign && \
touch alignment/FD25672365_2-2130752/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/FD25672365_2-2130752.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/FD25672365_2-2130752.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672365_2-2130752/FD25672365_2-2130752.sorted.realigned.bam
gatk_indel_realigner.FD25672365_2-2130752.899f2e0a1d7988810045196bab741e70.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_22_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_22_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_23_JOB_ID: gatk_indel_realigner.FD25672366_2-2130754
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672366_2-2130754
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672366_2-2130754.8d4872eb74df71524b1633fc148abf43.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672366_2-2130754.8d4872eb74df71524b1633fc148abf43.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672366_2-2130754/realign && \
touch alignment/FD25672366_2-2130754/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/FD25672366_2-2130754.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/FD25672366_2-2130754.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672366_2-2130754/FD25672366_2-2130754.sorted.realigned.bam
gatk_indel_realigner.FD25672366_2-2130754.8d4872eb74df71524b1633fc148abf43.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_23_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_23_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_24_JOB_ID: gatk_indel_realigner.FD25672367_2-2130756
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672367_2-2130756
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672367_2-2130756.4c94aad029d40ab2eb2ad7616b39ae3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672367_2-2130756.4c94aad029d40ab2eb2ad7616b39ae3b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672367_2-2130756/realign && \
touch alignment/FD25672367_2-2130756/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/FD25672367_2-2130756.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/FD25672367_2-2130756.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672367_2-2130756/FD25672367_2-2130756.sorted.realigned.bam
gatk_indel_realigner.FD25672367_2-2130756.4c94aad029d40ab2eb2ad7616b39ae3b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_24_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_24_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_25_JOB_ID: gatk_indel_realigner.FD25672368_2-2130758
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25672368_2-2130758
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25672368_2-2130758.a49abd6f5317a34d66642620ddc9e1e5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25672368_2-2130758.a49abd6f5317a34d66642620ddc9e1e5.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25672368_2-2130758/realign && \
touch alignment/FD25672368_2-2130758/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/FD25672368_2-2130758.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/FD25672368_2-2130758.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25672368_2-2130758/FD25672368_2-2130758.sorted.realigned.bam
gatk_indel_realigner.FD25672368_2-2130758.a49abd6f5317a34d66642620ddc9e1e5.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_25_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_25_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_26_JOB_ID: gatk_indel_realigner.FD25673014_2-2130766
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25673014_2-2130766
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25673014_2-2130766.0a006867871656177fb9f4f79ab52602.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25673014_2-2130766.0a006867871656177fb9f4f79ab52602.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25673014_2-2130766/realign && \
touch alignment/FD25673014_2-2130766/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/FD25673014_2-2130766.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/FD25673014_2-2130766.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25673014_2-2130766/FD25673014_2-2130766.sorted.realigned.bam
gatk_indel_realigner.FD25673014_2-2130766.0a006867871656177fb9f4f79ab52602.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_26_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_26_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_27_JOB_ID: gatk_indel_realigner.FD25699810_2-2130688
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699810_2-2130688
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699810_2-2130688.a3aa6e76d136d628f6e13497a4d9e1f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699810_2-2130688.a3aa6e76d136d628f6e13497a4d9e1f0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699810_2-2130688/realign && \
touch alignment/FD25699810_2-2130688/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/FD25699810_2-2130688.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/FD25699810_2-2130688.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699810_2-2130688/FD25699810_2-2130688.sorted.realigned.bam
gatk_indel_realigner.FD25699810_2-2130688.a3aa6e76d136d628f6e13497a4d9e1f0.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_27_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_27_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_28_JOB_ID: gatk_indel_realigner.FD25699821_2-2130573
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699821_2-2130573
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699821_2-2130573.d50b6b727cd89615f0f0d2626dfed856.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699821_2-2130573.d50b6b727cd89615f0f0d2626dfed856.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699821_2-2130573/realign && \
touch alignment/FD25699821_2-2130573/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/FD25699821_2-2130573.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/FD25699821_2-2130573.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699821_2-2130573/FD25699821_2-2130573.sorted.realigned.bam
gatk_indel_realigner.FD25699821_2-2130573.d50b6b727cd89615f0f0d2626dfed856.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_28_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_28_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_29_JOB_ID: gatk_indel_realigner.FD25699822_2-2130576
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699822_2-2130576
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699822_2-2130576.5da6fd946b007262c3b66c539e67edd7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699822_2-2130576.5da6fd946b007262c3b66c539e67edd7.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699822_2-2130576/realign && \
touch alignment/FD25699822_2-2130576/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/FD25699822_2-2130576.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/FD25699822_2-2130576.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699822_2-2130576/FD25699822_2-2130576.sorted.realigned.bam
gatk_indel_realigner.FD25699822_2-2130576.5da6fd946b007262c3b66c539e67edd7.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_29_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_29_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_30_JOB_ID: gatk_indel_realigner.FD25699823_2-2130579
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699823_2-2130579
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699823_2-2130579.d94bccb7dfca47cad28827e5764b187f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699823_2-2130579.d94bccb7dfca47cad28827e5764b187f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699823_2-2130579/realign && \
touch alignment/FD25699823_2-2130579/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/FD25699823_2-2130579.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/FD25699823_2-2130579.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699823_2-2130579/FD25699823_2-2130579.sorted.realigned.bam
gatk_indel_realigner.FD25699823_2-2130579.d94bccb7dfca47cad28827e5764b187f.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_30_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_30_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_31_JOB_ID: gatk_indel_realigner.FD25699829_2-2130591
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699829_2-2130591
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699829_2-2130591.df4575d6cdcde769712481e8f47f5732.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699829_2-2130591.df4575d6cdcde769712481e8f47f5732.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699829_2-2130591/realign && \
touch alignment/FD25699829_2-2130591/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/FD25699829_2-2130591.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/FD25699829_2-2130591.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699829_2-2130591/FD25699829_2-2130591.sorted.realigned.bam
gatk_indel_realigner.FD25699829_2-2130591.df4575d6cdcde769712481e8f47f5732.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_31_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_31_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_32_JOB_ID: gatk_indel_realigner.FD25699830_2-2130594
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699830_2-2130594
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699830_2-2130594.11f9b7d792348bd3e6c7ad469ece3da7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699830_2-2130594.11f9b7d792348bd3e6c7ad469ece3da7.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699830_2-2130594/realign && \
touch alignment/FD25699830_2-2130594/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/FD25699830_2-2130594.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/FD25699830_2-2130594.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699830_2-2130594/FD25699830_2-2130594.sorted.realigned.bam
gatk_indel_realigner.FD25699830_2-2130594.11f9b7d792348bd3e6c7ad469ece3da7.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_32_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_32_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_32_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_33_JOB_ID: gatk_indel_realigner.FD25699831_2-2130582
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699831_2-2130582
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699831_2-2130582.744881e4f4b46191ca80c9a9be013ab5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699831_2-2130582.744881e4f4b46191ca80c9a9be013ab5.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699831_2-2130582/realign && \
touch alignment/FD25699831_2-2130582/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/FD25699831_2-2130582.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/FD25699831_2-2130582.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699831_2-2130582/FD25699831_2-2130582.sorted.realigned.bam
gatk_indel_realigner.FD25699831_2-2130582.744881e4f4b46191ca80c9a9be013ab5.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_33_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_33_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_33_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_34_JOB_ID: gatk_indel_realigner.FD25699832_2-2130585
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699832_2-2130585
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699832_2-2130585.58db8df0b0f0d80f0638cf62e8574a0e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699832_2-2130585.58db8df0b0f0d80f0638cf62e8574a0e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699832_2-2130585/realign && \
touch alignment/FD25699832_2-2130585/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/FD25699832_2-2130585.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/FD25699832_2-2130585.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699832_2-2130585/FD25699832_2-2130585.sorted.realigned.bam
gatk_indel_realigner.FD25699832_2-2130585.58db8df0b0f0d80f0638cf62e8574a0e.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_34_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_34_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_34_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_35_JOB_ID: gatk_indel_realigner.FD25699833_2-2130588
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699833_2-2130588
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699833_2-2130588.c20d35063df5c20f0c26098115278343.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699833_2-2130588.c20d35063df5c20f0c26098115278343.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699833_2-2130588/realign && \
touch alignment/FD25699833_2-2130588/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/FD25699833_2-2130588.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/FD25699833_2-2130588.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699833_2-2130588/FD25699833_2-2130588.sorted.realigned.bam
gatk_indel_realigner.FD25699833_2-2130588.c20d35063df5c20f0c26098115278343.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_35_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_35_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_35_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_36_JOB_ID: gatk_indel_realigner.FD25699894_2-2130670
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699894_2-2130670
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699894_2-2130670.716743827cd917586a27f40cea3ebc53.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699894_2-2130670.716743827cd917586a27f40cea3ebc53.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699894_2-2130670/realign && \
touch alignment/FD25699894_2-2130670/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/FD25699894_2-2130670.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/FD25699894_2-2130670.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/FD25699894_2-2130670.sorted.realigned.bam
gatk_indel_realigner.FD25699894_2-2130670.716743827cd917586a27f40cea3ebc53.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_36_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_36_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_36_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_37_JOB_ID: gatk_indel_realigner.FD25699895_2-2130673
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699895_2-2130673
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699895_2-2130673.e1665219c7003309387cf0d24851d9c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699895_2-2130673.e1665219c7003309387cf0d24851d9c9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699895_2-2130673/realign && \
touch alignment/FD25699895_2-2130673/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/FD25699895_2-2130673.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/FD25699895_2-2130673.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699895_2-2130673/FD25699895_2-2130673.sorted.realigned.bam
gatk_indel_realigner.FD25699895_2-2130673.e1665219c7003309387cf0d24851d9c9.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_37_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_37_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_37_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_38_JOB_ID: gatk_indel_realigner.FD25699897_2-2130685
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699897_2-2130685
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699897_2-2130685.35059a86df0c8913f5afc1d57100a71a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699897_2-2130685.35059a86df0c8913f5afc1d57100a71a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699897_2-2130685/realign && \
touch alignment/FD25699897_2-2130685/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/FD25699897_2-2130685.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/FD25699897_2-2130685.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699897_2-2130685/FD25699897_2-2130685.sorted.realigned.bam
gatk_indel_realigner.FD25699897_2-2130685.35059a86df0c8913f5afc1d57100a71a.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_38_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_38_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_38_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_39_JOB_ID: gatk_indel_realigner.FD25699898_2-2130676
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699898_2-2130676
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699898_2-2130676.64a1151fad85f88cc86a32fca40350e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699898_2-2130676.64a1151fad85f88cc86a32fca40350e2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699898_2-2130676/realign && \
touch alignment/FD25699898_2-2130676/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/FD25699898_2-2130676.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/FD25699898_2-2130676.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699898_2-2130676/FD25699898_2-2130676.sorted.realigned.bam
gatk_indel_realigner.FD25699898_2-2130676.64a1151fad85f88cc86a32fca40350e2.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_39_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_39_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_39_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_40_JOB_ID: gatk_indel_realigner.FD25699899_2-2130679
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699899_2-2130679
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699899_2-2130679.1e650e997683a638c3cb8954b78aa24d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699899_2-2130679.1e650e997683a638c3cb8954b78aa24d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699899_2-2130679/realign && \
touch alignment/FD25699899_2-2130679/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/FD25699899_2-2130679.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/FD25699899_2-2130679.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699899_2-2130679/FD25699899_2-2130679.sorted.realigned.bam
gatk_indel_realigner.FD25699899_2-2130679.1e650e997683a638c3cb8954b78aa24d.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_40_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_40_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_40_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_41_JOB_ID: gatk_indel_realigner.FD25699900_2-2130682
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699900_2-2130682
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699900_2-2130682.6ffcc3678492ea9cf7b788d4aade73d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699900_2-2130682.6ffcc3678492ea9cf7b788d4aade73d6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699900_2-2130682/realign && \
touch alignment/FD25699900_2-2130682/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/FD25699900_2-2130682.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/FD25699900_2-2130682.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699900_2-2130682/FD25699900_2-2130682.sorted.realigned.bam
gatk_indel_realigner.FD25699900_2-2130682.6ffcc3678492ea9cf7b788d4aade73d6.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_41_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_41_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_41_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_42_JOB_ID: gatk_indel_realigner.FD25699901_2-2134255
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699901_2-2134255
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699901_2-2134255.edf0cb7bd1e3cbedc7beb13597696069.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699901_2-2134255.edf0cb7bd1e3cbedc7beb13597696069.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699901_2-2134255/realign && \
touch alignment/FD25699901_2-2134255/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/FD25699901_2-2134255.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/FD25699901_2-2134255.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699901_2-2134255/FD25699901_2-2134255.sorted.realigned.bam
gatk_indel_realigner.FD25699901_2-2134255.edf0cb7bd1e3cbedc7beb13597696069.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_42_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_42_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_42_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_43_JOB_ID: gatk_indel_realigner.FD25699902_2-2134257
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699902_2-2134257
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699902_2-2134257.987c66d5bfb3e2a490ea5d5b7848017e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699902_2-2134257.987c66d5bfb3e2a490ea5d5b7848017e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699902_2-2134257/realign && \
touch alignment/FD25699902_2-2134257/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/FD25699902_2-2134257.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/FD25699902_2-2134257.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699902_2-2134257/FD25699902_2-2134257.sorted.realigned.bam
gatk_indel_realigner.FD25699902_2-2134257.987c66d5bfb3e2a490ea5d5b7848017e.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_43_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_43_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_43_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_44_JOB_ID: gatk_indel_realigner.FD25699903_2-2134259
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699903_2-2134259
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699903_2-2134259.1fc45e3569590703452b00c13b3d263e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699903_2-2134259.1fc45e3569590703452b00c13b3d263e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699903_2-2134259/realign && \
touch alignment/FD25699903_2-2134259/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/FD25699903_2-2134259.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/FD25699903_2-2134259.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699903_2-2134259/FD25699903_2-2134259.sorted.realigned.bam
gatk_indel_realigner.FD25699903_2-2134259.1fc45e3569590703452b00c13b3d263e.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_44_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_44_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_44_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_45_JOB_ID: gatk_indel_realigner.FD25699904_2-2134261
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699904_2-2134261
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699904_2-2134261.24f594c5c1eb10f71116b058ab5b7032.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699904_2-2134261.24f594c5c1eb10f71116b058ab5b7032.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699904_2-2134261/realign && \
touch alignment/FD25699904_2-2134261/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/FD25699904_2-2134261.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/FD25699904_2-2134261.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699904_2-2134261/FD25699904_2-2134261.sorted.realigned.bam
gatk_indel_realigner.FD25699904_2-2134261.24f594c5c1eb10f71116b058ab5b7032.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_45_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_45_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_45_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_46_JOB_ID: gatk_indel_realigner.FD25699905_2-2134099
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699905_2-2134099
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699905_2-2134099.6a7cef0c63a692a26ab2964b5660966d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699905_2-2134099.6a7cef0c63a692a26ab2964b5660966d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699905_2-2134099/realign && \
touch alignment/FD25699905_2-2134099/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/FD25699905_2-2134099.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/FD25699905_2-2134099.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699905_2-2134099/FD25699905_2-2134099.sorted.realigned.bam
gatk_indel_realigner.FD25699905_2-2134099.6a7cef0c63a692a26ab2964b5660966d.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_46_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_46_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_46_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_47_JOB_ID: gatk_indel_realigner.FD25699906_2-2134251
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699906_2-2134251
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699906_2-2134251.b24a389a412026caf91dbe0f9298ba3e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699906_2-2134251.b24a389a412026caf91dbe0f9298ba3e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699906_2-2134251/realign && \
touch alignment/FD25699906_2-2134251/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/FD25699906_2-2134251.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/FD25699906_2-2134251.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699906_2-2134251/FD25699906_2-2134251.sorted.realigned.bam
gatk_indel_realigner.FD25699906_2-2134251.b24a389a412026caf91dbe0f9298ba3e.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_47_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_47_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_47_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_48_JOB_ID: gatk_indel_realigner.FD25699907_2-2134253
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699907_2-2134253
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699907_2-2134253.a8acda5a6d6bca47e3f5940566284d73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699907_2-2134253.a8acda5a6d6bca47e3f5940566284d73.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699907_2-2134253/realign && \
touch alignment/FD25699907_2-2134253/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/FD25699907_2-2134253.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/FD25699907_2-2134253.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699907_2-2134253/FD25699907_2-2134253.sorted.realigned.bam
gatk_indel_realigner.FD25699907_2-2134253.a8acda5a6d6bca47e3f5940566284d73.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_48_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_48_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_48_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_49_JOB_ID: gatk_indel_realigner.FD25699909_2-2130590
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699909_2-2130590
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699909_2-2130590.dbd73af84db0d415bb028368f867dfb6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699909_2-2130590.dbd73af84db0d415bb028368f867dfb6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699909_2-2130590/realign && \
touch alignment/FD25699909_2-2130590/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/FD25699909_2-2130590.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/FD25699909_2-2130590.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699909_2-2130590/FD25699909_2-2130590.sorted.realigned.bam
gatk_indel_realigner.FD25699909_2-2130590.dbd73af84db0d415bb028368f867dfb6.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_49_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_49_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_49_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_50_JOB_ID: gatk_indel_realigner.FD25699910_2-2130581
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699910_2-2130581
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699910_2-2130581.5e0a0d9a81a7c42ab856d6ab2d2badac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699910_2-2130581.5e0a0d9a81a7c42ab856d6ab2d2badac.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699910_2-2130581/realign && \
touch alignment/FD25699910_2-2130581/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/FD25699910_2-2130581.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/FD25699910_2-2130581.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699910_2-2130581/FD25699910_2-2130581.sorted.realigned.bam
gatk_indel_realigner.FD25699910_2-2130581.5e0a0d9a81a7c42ab856d6ab2d2badac.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_50_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_50_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_50_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_51_JOB_ID: gatk_indel_realigner.FD25699911_2-2130584
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699911_2-2130584
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699911_2-2130584.c0b80a9d2b5876ceb794353fb4c1534e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699911_2-2130584.c0b80a9d2b5876ceb794353fb4c1534e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699911_2-2130584/realign && \
touch alignment/FD25699911_2-2130584/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/FD25699911_2-2130584.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/FD25699911_2-2130584.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699911_2-2130584/FD25699911_2-2130584.sorted.realigned.bam
gatk_indel_realigner.FD25699911_2-2130584.c0b80a9d2b5876ceb794353fb4c1534e.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_51_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_51_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_51_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_52_JOB_ID: gatk_indel_realigner.FD25699912_2-2130587
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699912_2-2130587
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699912_2-2130587.c4b70e18f4ce08e32c273192df94d9c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699912_2-2130587.c4b70e18f4ce08e32c273192df94d9c3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699912_2-2130587/realign && \
touch alignment/FD25699912_2-2130587/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/FD25699912_2-2130587.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/FD25699912_2-2130587.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699912_2-2130587/FD25699912_2-2130587.sorted.realigned.bam
gatk_indel_realigner.FD25699912_2-2130587.c4b70e18f4ce08e32c273192df94d9c3.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_52_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_52_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_52_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_53_JOB_ID: gatk_indel_realigner.FD25699913_2-2130578
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699913_2-2130578
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699913_2-2130578.e2b820d7c1cc2f81339158308143561d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699913_2-2130578.e2b820d7c1cc2f81339158308143561d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699913_2-2130578/realign && \
touch alignment/FD25699913_2-2130578/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/FD25699913_2-2130578.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/FD25699913_2-2130578.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699913_2-2130578/FD25699913_2-2130578.sorted.realigned.bam
gatk_indel_realigner.FD25699913_2-2130578.e2b820d7c1cc2f81339158308143561d.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_53_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_53_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_53_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_54_JOB_ID: gatk_indel_realigner.FD25699914_2-2130572
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699914_2-2130572
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699914_2-2130572.9cea37090f25b1c8da2f8ac54ef95526.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699914_2-2130572.9cea37090f25b1c8da2f8ac54ef95526.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699914_2-2130572/realign && \
touch alignment/FD25699914_2-2130572/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/FD25699914_2-2130572.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/FD25699914_2-2130572.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699914_2-2130572/FD25699914_2-2130572.sorted.realigned.bam
gatk_indel_realigner.FD25699914_2-2130572.9cea37090f25b1c8da2f8ac54ef95526.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_54_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_54_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_54_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_55_JOB_ID: gatk_indel_realigner.FD25699915_2-2130575
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699915_2-2130575
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699915_2-2130575.0cd41665dcad3921b50979b9fa664f85.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699915_2-2130575.0cd41665dcad3921b50979b9fa664f85.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699915_2-2130575/realign && \
touch alignment/FD25699915_2-2130575/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/FD25699915_2-2130575.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/FD25699915_2-2130575.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699915_2-2130575/FD25699915_2-2130575.sorted.realigned.bam
gatk_indel_realigner.FD25699915_2-2130575.0cd41665dcad3921b50979b9fa664f85.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_55_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_55_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_55_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_56_JOB_ID: gatk_indel_realigner.FD25699916_2-2130669
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699916_2-2130669
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699916_2-2130669.656300ec35b27b19ca0872668e132c36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699916_2-2130669.656300ec35b27b19ca0872668e132c36.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699916_2-2130669/realign && \
touch alignment/FD25699916_2-2130669/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/FD25699916_2-2130669.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/FD25699916_2-2130669.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699916_2-2130669/FD25699916_2-2130669.sorted.realigned.bam
gatk_indel_realigner.FD25699916_2-2130669.656300ec35b27b19ca0872668e132c36.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_56_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_56_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_56_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_57_JOB_ID: gatk_indel_realigner.FD25699917_2-2130663
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699917_2-2130663
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699917_2-2130663.6bd956db7726cadbcec5a7b5618c0dbd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699917_2-2130663.6bd956db7726cadbcec5a7b5618c0dbd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699917_2-2130663/realign && \
touch alignment/FD25699917_2-2130663/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/FD25699917_2-2130663.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/FD25699917_2-2130663.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699917_2-2130663/FD25699917_2-2130663.sorted.realigned.bam
gatk_indel_realigner.FD25699917_2-2130663.6bd956db7726cadbcec5a7b5618c0dbd.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_57_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_57_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_57_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_58_JOB_ID: gatk_indel_realigner.FD25699918_2-2130666
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699918_2-2130666
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699918_2-2130666.1ef4114c4a6c7014c468d0c16e8cafd8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699918_2-2130666.1ef4114c4a6c7014c468d0c16e8cafd8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699918_2-2130666/realign && \
touch alignment/FD25699918_2-2130666/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/FD25699918_2-2130666.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/FD25699918_2-2130666.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699918_2-2130666/FD25699918_2-2130666.sorted.realigned.bam
gatk_indel_realigner.FD25699918_2-2130666.1ef4114c4a6c7014c468d0c16e8cafd8.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_58_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_58_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_58_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_59_JOB_ID: gatk_indel_realigner.FD25699919_2-2130657
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699919_2-2130657
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699919_2-2130657.b93ee075a72cb46dd11929699efca6f4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699919_2-2130657.b93ee075a72cb46dd11929699efca6f4.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699919_2-2130657/realign && \
touch alignment/FD25699919_2-2130657/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/FD25699919_2-2130657.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/FD25699919_2-2130657.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699919_2-2130657/FD25699919_2-2130657.sorted.realigned.bam
gatk_indel_realigner.FD25699919_2-2130657.b93ee075a72cb46dd11929699efca6f4.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_59_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_59_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_59_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_60_JOB_ID: gatk_indel_realigner.FD25699920_2-2130660
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699920_2-2130660
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699920_2-2130660.91beef75feca20150c1a611684a77f08.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699920_2-2130660.91beef75feca20150c1a611684a77f08.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699920_2-2130660/realign && \
touch alignment/FD25699920_2-2130660/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/FD25699920_2-2130660.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/FD25699920_2-2130660.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699920_2-2130660/FD25699920_2-2130660.sorted.realigned.bam
gatk_indel_realigner.FD25699920_2-2130660.91beef75feca20150c1a611684a77f08.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_60_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_60_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_60_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_61_JOB_ID: gatk_indel_realigner.FD25699921_2-2130598
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699921_2-2130598
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699921_2-2130598.c0a62347f18da613cac06c0d257a5d98.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699921_2-2130598.c0a62347f18da613cac06c0d257a5d98.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699921_2-2130598/realign && \
touch alignment/FD25699921_2-2130598/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/FD25699921_2-2130598.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/FD25699921_2-2130598.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699921_2-2130598/FD25699921_2-2130598.sorted.realigned.bam
gatk_indel_realigner.FD25699921_2-2130598.c0a62347f18da613cac06c0d257a5d98.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_61_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_61_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_61_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_62_JOB_ID: gatk_indel_realigner.FD25699922_2-2130651
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699922_2-2130651
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699922_2-2130651.0c8a47be9e84285cae0a9fab24cfc88f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699922_2-2130651.0c8a47be9e84285cae0a9fab24cfc88f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699922_2-2130651/realign && \
touch alignment/FD25699922_2-2130651/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/FD25699922_2-2130651.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/FD25699922_2-2130651.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699922_2-2130651/FD25699922_2-2130651.sorted.realigned.bam
gatk_indel_realigner.FD25699922_2-2130651.0c8a47be9e84285cae0a9fab24cfc88f.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_62_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_62_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_62_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_63_JOB_ID: gatk_indel_realigner.FD25699923_2-2130654
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699923_2-2130654
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699923_2-2130654.0d2d72c7301d49781d6cb159a4988082.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699923_2-2130654.0d2d72c7301d49781d6cb159a4988082.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699923_2-2130654/realign && \
touch alignment/FD25699923_2-2130654/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/FD25699923_2-2130654.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/FD25699923_2-2130654.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699923_2-2130654/FD25699923_2-2130654.sorted.realigned.bam
gatk_indel_realigner.FD25699923_2-2130654.0d2d72c7301d49781d6cb159a4988082.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_63_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_63_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_63_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_64_JOB_ID: gatk_indel_realigner.FD25699941_2-2130674
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699941_2-2130674
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699941_2-2130674.ab37faf35fd4c4348b51a3350fd63686.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699941_2-2130674.ab37faf35fd4c4348b51a3350fd63686.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699941_2-2130674/realign && \
touch alignment/FD25699941_2-2130674/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/FD25699941_2-2130674.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/FD25699941_2-2130674.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699941_2-2130674/FD25699941_2-2130674.sorted.realigned.bam
gatk_indel_realigner.FD25699941_2-2130674.ab37faf35fd4c4348b51a3350fd63686.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_64_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_64_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_64_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_65_JOB_ID: gatk_indel_realigner.FD25699942_2-2130677
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699942_2-2130677
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699942_2-2130677.e051c05af00b31b588e1608080926d0b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699942_2-2130677.e051c05af00b31b588e1608080926d0b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699942_2-2130677/realign && \
touch alignment/FD25699942_2-2130677/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/FD25699942_2-2130677.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/FD25699942_2-2130677.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699942_2-2130677/FD25699942_2-2130677.sorted.realigned.bam
gatk_indel_realigner.FD25699942_2-2130677.e051c05af00b31b588e1608080926d0b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_65_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_65_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_65_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_66_JOB_ID: gatk_indel_realigner.FD25699952_2-2130689
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699952_2-2130689
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699952_2-2130689.6dfb24f1453d311ddc86a3c6e583d617.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699952_2-2130689.6dfb24f1453d311ddc86a3c6e583d617.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699952_2-2130689/realign && \
touch alignment/FD25699952_2-2130689/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/FD25699952_2-2130689.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/FD25699952_2-2130689.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699952_2-2130689/FD25699952_2-2130689.sorted.realigned.bam
gatk_indel_realigner.FD25699952_2-2130689.6dfb24f1453d311ddc86a3c6e583d617.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_66_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_66_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_66_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_67_JOB_ID: gatk_indel_realigner.FD25699953_2-2130680
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699953_2-2130680
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699953_2-2130680.4c1a90ced455760d880e21f5912a434a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699953_2-2130680.4c1a90ced455760d880e21f5912a434a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699953_2-2130680/realign && \
touch alignment/FD25699953_2-2130680/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/FD25699953_2-2130680.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/FD25699953_2-2130680.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699953_2-2130680/FD25699953_2-2130680.sorted.realigned.bam
gatk_indel_realigner.FD25699953_2-2130680.4c1a90ced455760d880e21f5912a434a.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_67_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_67_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_67_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_68_JOB_ID: gatk_indel_realigner.FD25699954_2-2130683
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699954_2-2130683
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699954_2-2130683.db19cc507b4f95072106c41186a13544.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699954_2-2130683.db19cc507b4f95072106c41186a13544.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699954_2-2130683/realign && \
touch alignment/FD25699954_2-2130683/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/FD25699954_2-2130683.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/FD25699954_2-2130683.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699954_2-2130683/FD25699954_2-2130683.sorted.realigned.bam
gatk_indel_realigner.FD25699954_2-2130683.db19cc507b4f95072106c41186a13544.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_68_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_68_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_68_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_69_JOB_ID: gatk_indel_realigner.FD25699955_2-2130686
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699955_2-2130686
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699955_2-2130686.c5a81882c1089017673cc97775bb8bc3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699955_2-2130686.c5a81882c1089017673cc97775bb8bc3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699955_2-2130686/realign && \
touch alignment/FD25699955_2-2130686/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/FD25699955_2-2130686.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/FD25699955_2-2130686.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699955_2-2130686/FD25699955_2-2130686.sorted.realigned.bam
gatk_indel_realigner.FD25699955_2-2130686.c5a81882c1089017673cc97775bb8bc3.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_69_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_69_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_69_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_70_JOB_ID: gatk_indel_realigner.FD25699956_2-2130671
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699956_2-2130671
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699956_2-2130671.ce299994111ad1ad762c386e78c66498.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699956_2-2130671.ce299994111ad1ad762c386e78c66498.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699956_2-2130671/realign && \
touch alignment/FD25699956_2-2130671/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/FD25699956_2-2130671.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/FD25699956_2-2130671.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699956_2-2130671/FD25699956_2-2130671.sorted.realigned.bam
gatk_indel_realigner.FD25699956_2-2130671.ce299994111ad1ad762c386e78c66498.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_70_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_70_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_70_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_71_JOB_ID: gatk_indel_realigner.FD25699964_2-2130692
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699964_2-2130692
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699964_2-2130692.5aad56ffc71003c2d82f991343a6e68b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699964_2-2130692.5aad56ffc71003c2d82f991343a6e68b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699964_2-2130692/realign && \
touch alignment/FD25699964_2-2130692/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/FD25699964_2-2130692.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/FD25699964_2-2130692.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699964_2-2130692/FD25699964_2-2130692.sorted.realigned.bam
gatk_indel_realigner.FD25699964_2-2130692.5aad56ffc71003c2d82f991343a6e68b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_71_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_71_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_71_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_72_JOB_ID: gatk_indel_realigner.FD25699986_2-2130687
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699986_2-2130687
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699986_2-2130687.1263b7f9ace5b40df640e60e81f85349.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699986_2-2130687.1263b7f9ace5b40df640e60e81f85349.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699986_2-2130687/realign && \
touch alignment/FD25699986_2-2130687/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/FD25699986_2-2130687.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/FD25699986_2-2130687.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699986_2-2130687/FD25699986_2-2130687.sorted.realigned.bam
gatk_indel_realigner.FD25699986_2-2130687.1263b7f9ace5b40df640e60e81f85349.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_72_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_72_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_72_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_73_JOB_ID: gatk_indel_realigner.FD25699987_2-2130690
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699987_2-2130690
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699987_2-2130690.e5c043f4335371a4c108e70ff5cad97c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699987_2-2130690.e5c043f4335371a4c108e70ff5cad97c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699987_2-2130690/realign && \
touch alignment/FD25699987_2-2130690/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/FD25699987_2-2130690.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/FD25699987_2-2130690.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699987_2-2130690/FD25699987_2-2130690.sorted.realigned.bam
gatk_indel_realigner.FD25699987_2-2130690.e5c043f4335371a4c108e70ff5cad97c.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_73_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_73_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_73_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_74_JOB_ID: gatk_indel_realigner.FD25699989_2-2130672
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699989_2-2130672
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699989_2-2130672.a7df7f3463b934c54121d890f43ebd7c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699989_2-2130672.a7df7f3463b934c54121d890f43ebd7c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699989_2-2130672/realign && \
touch alignment/FD25699989_2-2130672/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/FD25699989_2-2130672.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/FD25699989_2-2130672.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699989_2-2130672/FD25699989_2-2130672.sorted.realigned.bam
gatk_indel_realigner.FD25699989_2-2130672.a7df7f3463b934c54121d890f43ebd7c.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_74_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_74_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_74_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_75_JOB_ID: gatk_indel_realigner.FD25699990_2-2130675
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699990_2-2130675
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699990_2-2130675.1d14091f89d7e5ae17081c28706273aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699990_2-2130675.1d14091f89d7e5ae17081c28706273aa.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699990_2-2130675/realign && \
touch alignment/FD25699990_2-2130675/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/FD25699990_2-2130675.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/FD25699990_2-2130675.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699990_2-2130675/FD25699990_2-2130675.sorted.realigned.bam
gatk_indel_realigner.FD25699990_2-2130675.1d14091f89d7e5ae17081c28706273aa.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_75_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_75_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_75_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_76_JOB_ID: gatk_indel_realigner.FD25699991_2-2130678
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699991_2-2130678
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699991_2-2130678.0355f732be6c815684f392bd9103d4fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699991_2-2130678.0355f732be6c815684f392bd9103d4fc.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699991_2-2130678/realign && \
touch alignment/FD25699991_2-2130678/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/FD25699991_2-2130678.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/FD25699991_2-2130678.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699991_2-2130678/FD25699991_2-2130678.sorted.realigned.bam
gatk_indel_realigner.FD25699991_2-2130678.0355f732be6c815684f392bd9103d4fc.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_76_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_76_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_76_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_77_JOB_ID: gatk_indel_realigner.FD25699992_2-2130681
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699992_2-2130681
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699992_2-2130681.cef9b5fa1e30616673b1684de3e0a039.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699992_2-2130681.cef9b5fa1e30616673b1684de3e0a039.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699992_2-2130681/realign && \
touch alignment/FD25699992_2-2130681/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/FD25699992_2-2130681.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/FD25699992_2-2130681.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699992_2-2130681/FD25699992_2-2130681.sorted.realigned.bam
gatk_indel_realigner.FD25699992_2-2130681.cef9b5fa1e30616673b1684de3e0a039.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_77_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_77_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_77_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_78_JOB_ID: gatk_indel_realigner.FD25699993_2-2130684
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25699993_2-2130684
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25699993_2-2130684.85272a687a5617715450c1236ba517a0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25699993_2-2130684.85272a687a5617715450c1236ba517a0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25699993_2-2130684/realign && \
touch alignment/FD25699993_2-2130684/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/FD25699993_2-2130684.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/FD25699993_2-2130684.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699993_2-2130684/FD25699993_2-2130684.sorted.realigned.bam
gatk_indel_realigner.FD25699993_2-2130684.85272a687a5617715450c1236ba517a0.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_78_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_78_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_78_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_79_JOB_ID: gatk_indel_realigner.FD25700005_2-2130662
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700005_2-2130662
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700005_2-2130662.e18765b572cb9c255b2d91fdb9a8361b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700005_2-2130662.e18765b572cb9c255b2d91fdb9a8361b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700005_2-2130662/realign && \
touch alignment/FD25700005_2-2130662/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/FD25700005_2-2130662.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/FD25700005_2-2130662.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700005_2-2130662/FD25700005_2-2130662.sorted.realigned.bam
gatk_indel_realigner.FD25700005_2-2130662.e18765b572cb9c255b2d91fdb9a8361b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_79_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_79_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_79_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_80_JOB_ID: gatk_indel_realigner.FD25700006_2-2130665
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700006_2-2130665
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700006_2-2130665.8780ed2a5b3ee222a87a2105af248268.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700006_2-2130665.8780ed2a5b3ee222a87a2105af248268.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700006_2-2130665/realign && \
touch alignment/FD25700006_2-2130665/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/FD25700006_2-2130665.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/FD25700006_2-2130665.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700006_2-2130665/FD25700006_2-2130665.sorted.realigned.bam
gatk_indel_realigner.FD25700006_2-2130665.8780ed2a5b3ee222a87a2105af248268.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_80_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_80_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_80_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_81_JOB_ID: gatk_indel_realigner.FD25700007_2-2130668
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700007_2-2130668
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700007_2-2130668.b03a133013b3e14d0d8c0c29a4c4521b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700007_2-2130668.b03a133013b3e14d0d8c0c29a4c4521b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700007_2-2130668/realign && \
touch alignment/FD25700007_2-2130668/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/FD25700007_2-2130668.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/FD25700007_2-2130668.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700007_2-2130668/FD25700007_2-2130668.sorted.realigned.bam
gatk_indel_realigner.FD25700007_2-2130668.b03a133013b3e14d0d8c0c29a4c4521b.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_81_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_81_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_81_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_82_JOB_ID: gatk_indel_realigner.FD25700013_2-2134276
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700013_2-2134276
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700013_2-2134276.cf4e1745a9b8a5867d07583c267ff778.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700013_2-2134276.cf4e1745a9b8a5867d07583c267ff778.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700013_2-2134276/realign && \
touch alignment/FD25700013_2-2134276/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/FD25700013_2-2134276.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/FD25700013_2-2134276.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700013_2-2134276/FD25700013_2-2134276.sorted.realigned.bam
gatk_indel_realigner.FD25700013_2-2134276.cf4e1745a9b8a5867d07583c267ff778.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_82_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_82_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_82_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_83_JOB_ID: gatk_indel_realigner.FD25700014_2-2134278
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700014_2-2134278
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700014_2-2134278.395a15f947a32b5b6af272114219c42c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700014_2-2134278.395a15f947a32b5b6af272114219c42c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700014_2-2134278/realign && \
touch alignment/FD25700014_2-2134278/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/FD25700014_2-2134278.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/FD25700014_2-2134278.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700014_2-2134278/FD25700014_2-2134278.sorted.realigned.bam
gatk_indel_realigner.FD25700014_2-2134278.395a15f947a32b5b6af272114219c42c.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_83_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_83_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_83_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_84_JOB_ID: gatk_indel_realigner.FD25700015_2-2134280
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700015_2-2134280
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700015_2-2134280.31725a16347d3cf30be5aadcf445e435.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700015_2-2134280.31725a16347d3cf30be5aadcf445e435.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700015_2-2134280/realign && \
touch alignment/FD25700015_2-2134280/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/FD25700015_2-2134280.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/FD25700015_2-2134280.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700015_2-2134280/FD25700015_2-2134280.sorted.realigned.bam
gatk_indel_realigner.FD25700015_2-2134280.31725a16347d3cf30be5aadcf445e435.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_84_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_84_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_84_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_85_JOB_ID: gatk_indel_realigner.FD25700016_2-2130597
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700016_2-2130597
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700016_2-2130597.f0068cf060b86a8acdb91ab1d9419ee6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700016_2-2130597.f0068cf060b86a8acdb91ab1d9419ee6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700016_2-2130597/realign && \
touch alignment/FD25700016_2-2130597/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/FD25700016_2-2130597.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/FD25700016_2-2130597.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700016_2-2130597/FD25700016_2-2130597.sorted.realigned.bam
gatk_indel_realigner.FD25700016_2-2130597.f0068cf060b86a8acdb91ab1d9419ee6.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_85_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_85_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_85_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_86_JOB_ID: gatk_indel_realigner.FD25700017_2-2130600
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700017_2-2130600
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700017_2-2130600.b0edaa92f68565d365d91b127f59f0cd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700017_2-2130600.b0edaa92f68565d365d91b127f59f0cd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700017_2-2130600/realign && \
touch alignment/FD25700017_2-2130600/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/FD25700017_2-2130600.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/FD25700017_2-2130600.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700017_2-2130600/FD25700017_2-2130600.sorted.realigned.bam
gatk_indel_realigner.FD25700017_2-2130600.b0edaa92f68565d365d91b127f59f0cd.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_86_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_86_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_86_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_87_JOB_ID: gatk_indel_realigner.FD25700018_2-2130653
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700018_2-2130653
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700018_2-2130653.2378de5dc9266396d4cb31654f048253.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700018_2-2130653.2378de5dc9266396d4cb31654f048253.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700018_2-2130653/realign && \
touch alignment/FD25700018_2-2130653/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/FD25700018_2-2130653.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/FD25700018_2-2130653.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700018_2-2130653/FD25700018_2-2130653.sorted.realigned.bam
gatk_indel_realigner.FD25700018_2-2130653.2378de5dc9266396d4cb31654f048253.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_87_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_87_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_87_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_88_JOB_ID: gatk_indel_realigner.FD25700019_2-2130656
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700019_2-2130656
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700019_2-2130656.1cea362b4b6e932d415f1618bfd910a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700019_2-2130656.1cea362b4b6e932d415f1618bfd910a9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700019_2-2130656/realign && \
touch alignment/FD25700019_2-2130656/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/FD25700019_2-2130656.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/FD25700019_2-2130656.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700019_2-2130656/FD25700019_2-2130656.sorted.realigned.bam
gatk_indel_realigner.FD25700019_2-2130656.1cea362b4b6e932d415f1618bfd910a9.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_88_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_88_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_88_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_89_JOB_ID: gatk_indel_realigner.FD25700020_2-2130659
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700020_2-2130659
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700020_2-2130659.9a706b376197dbf8a0c467f3bd78da88.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700020_2-2130659.9a706b376197dbf8a0c467f3bd78da88.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700020_2-2130659/realign && \
touch alignment/FD25700020_2-2130659/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/FD25700020_2-2130659.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/FD25700020_2-2130659.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700020_2-2130659/FD25700020_2-2130659.sorted.realigned.bam
gatk_indel_realigner.FD25700020_2-2130659.9a706b376197dbf8a0c467f3bd78da88.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_89_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_89_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_89_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_90_JOB_ID: gatk_indel_realigner.FD25700042_2-2130574
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700042_2-2130574
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700042_2-2130574.3de4b8723d4893bbf3026c4a5f384d89.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700042_2-2130574.3de4b8723d4893bbf3026c4a5f384d89.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700042_2-2130574/realign && \
touch alignment/FD25700042_2-2130574/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/FD25700042_2-2130574.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/FD25700042_2-2130574.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700042_2-2130574/FD25700042_2-2130574.sorted.realigned.bam
gatk_indel_realigner.FD25700042_2-2130574.3de4b8723d4893bbf3026c4a5f384d89.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_90_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_90_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_90_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_91_JOB_ID: gatk_indel_realigner.FD25700043_2-2130577
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700043_2-2130577
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700043_2-2130577.6c56b5cde9334f7bb806066f9517f024.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700043_2-2130577.6c56b5cde9334f7bb806066f9517f024.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700043_2-2130577/realign && \
touch alignment/FD25700043_2-2130577/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/FD25700043_2-2130577.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/FD25700043_2-2130577.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700043_2-2130577/FD25700043_2-2130577.sorted.realigned.bam
gatk_indel_realigner.FD25700043_2-2130577.6c56b5cde9334f7bb806066f9517f024.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_91_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_91_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_91_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_92_JOB_ID: gatk_indel_realigner.FD25700044_2-2130580
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700044_2-2130580
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700044_2-2130580.1b3520d3c17e63d8c7c17e65f6038208.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700044_2-2130580.1b3520d3c17e63d8c7c17e65f6038208.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700044_2-2130580/realign && \
touch alignment/FD25700044_2-2130580/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/FD25700044_2-2130580.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/FD25700044_2-2130580.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700044_2-2130580/FD25700044_2-2130580.sorted.realigned.bam
gatk_indel_realigner.FD25700044_2-2130580.1b3520d3c17e63d8c7c17e65f6038208.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_92_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_92_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_92_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_93_JOB_ID: gatk_indel_realigner.FD25700045_2-2130592
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700045_2-2130592
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700045_2-2130592.3b0c5982a0915de5af7e549c260fcc28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700045_2-2130592.3b0c5982a0915de5af7e549c260fcc28.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700045_2-2130592/realign && \
touch alignment/FD25700045_2-2130592/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/FD25700045_2-2130592.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/FD25700045_2-2130592.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700045_2-2130592/FD25700045_2-2130592.sorted.realigned.bam
gatk_indel_realigner.FD25700045_2-2130592.3b0c5982a0915de5af7e549c260fcc28.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_93_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_93_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_93_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_94_JOB_ID: gatk_indel_realigner.FD25700046_2-2130595
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700046_2-2130595
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700046_2-2130595.5b9e80c775c2c4e5d12770c5e88f6f90.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700046_2-2130595.5b9e80c775c2c4e5d12770c5e88f6f90.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700046_2-2130595/realign && \
touch alignment/FD25700046_2-2130595/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/FD25700046_2-2130595.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/FD25700046_2-2130595.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700046_2-2130595/FD25700046_2-2130595.sorted.realigned.bam
gatk_indel_realigner.FD25700046_2-2130595.5b9e80c775c2c4e5d12770c5e88f6f90.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_94_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_94_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_94_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_95_JOB_ID: gatk_indel_realigner.FD25700047_2-2130583
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700047_2-2130583
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700047_2-2130583.17f53946d73c6ae9d2a1e04b40507b12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700047_2-2130583.17f53946d73c6ae9d2a1e04b40507b12.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700047_2-2130583/realign && \
touch alignment/FD25700047_2-2130583/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/FD25700047_2-2130583.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/FD25700047_2-2130583.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700047_2-2130583/FD25700047_2-2130583.sorted.realigned.bam
gatk_indel_realigner.FD25700047_2-2130583.17f53946d73c6ae9d2a1e04b40507b12.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_95_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_95_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_95_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_96_JOB_ID: gatk_indel_realigner.FD25700048_2-2130586
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700048_2-2130586
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700048_2-2130586.d68893994fecb937daa0f1434e86aaad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700048_2-2130586.d68893994fecb937daa0f1434e86aaad.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700048_2-2130586/realign && \
touch alignment/FD25700048_2-2130586/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/FD25700048_2-2130586.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/FD25700048_2-2130586.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700048_2-2130586/FD25700048_2-2130586.sorted.realigned.bam
gatk_indel_realigner.FD25700048_2-2130586.d68893994fecb937daa0f1434e86aaad.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_96_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_96_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_96_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_97_JOB_ID: gatk_indel_realigner.FD25700049_2-2130589
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700049_2-2130589
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700049_2-2130589.9bed17048ed3efb126a6f565e33c5942.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700049_2-2130589.9bed17048ed3efb126a6f565e33c5942.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700049_2-2130589/realign && \
touch alignment/FD25700049_2-2130589/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/FD25700049_2-2130589.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/FD25700049_2-2130589.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700049_2-2130589/FD25700049_2-2130589.sorted.realigned.bam
gatk_indel_realigner.FD25700049_2-2130589.9bed17048ed3efb126a6f565e33c5942.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_97_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_97_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_97_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_98_JOB_ID: gatk_indel_realigner.FD25700087_2-2134274
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700087_2-2134274
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700087_2-2134274.c2757a68a2a1bd443e1c819ed50b8055.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700087_2-2134274.c2757a68a2a1bd443e1c819ed50b8055.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700087_2-2134274/realign && \
touch alignment/FD25700087_2-2134274/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/FD25700087_2-2134274.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/FD25700087_2-2134274.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700087_2-2134274/FD25700087_2-2134274.sorted.realigned.bam
gatk_indel_realigner.FD25700087_2-2134274.c2757a68a2a1bd443e1c819ed50b8055.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_98_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_98_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_98_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_99_JOB_ID: gatk_indel_realigner.FD25700088_2-2134266
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700088_2-2134266
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700088_2-2134266.732f285eda7e66af6716236cafffa042.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700088_2-2134266.732f285eda7e66af6716236cafffa042.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700088_2-2134266/realign && \
touch alignment/FD25700088_2-2134266/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/FD25700088_2-2134266.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/FD25700088_2-2134266.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700088_2-2134266/FD25700088_2-2134266.sorted.realigned.bam
gatk_indel_realigner.FD25700088_2-2134266.732f285eda7e66af6716236cafffa042.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_99_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_99_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_99_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_100_JOB_ID: gatk_indel_realigner.FD25700089_2-2134268
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700089_2-2134268
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700089_2-2134268.358174190b7d8acd07a834424f93034a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700089_2-2134268.358174190b7d8acd07a834424f93034a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700089_2-2134268/realign && \
touch alignment/FD25700089_2-2134268/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/FD25700089_2-2134268.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/FD25700089_2-2134268.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700089_2-2134268/FD25700089_2-2134268.sorted.realigned.bam
gatk_indel_realigner.FD25700089_2-2134268.358174190b7d8acd07a834424f93034a.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_100_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_100_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_100_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_101_JOB_ID: gatk_indel_realigner.FD25700090_2-2134270
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700090_2-2134270
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700090_2-2134270.e6a36194382c43b855f25783a7fb29b8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700090_2-2134270.e6a36194382c43b855f25783a7fb29b8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700090_2-2134270/realign && \
touch alignment/FD25700090_2-2134270/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/FD25700090_2-2134270.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/FD25700090_2-2134270.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700090_2-2134270/FD25700090_2-2134270.sorted.realigned.bam
gatk_indel_realigner.FD25700090_2-2134270.e6a36194382c43b855f25783a7fb29b8.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_101_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_101_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_101_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_102_JOB_ID: gatk_indel_realigner.FD25700091_2-2134272
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700091_2-2134272
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700091_2-2134272.9ad49eea9a306256c2a3459d90ffc362.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700091_2-2134272.9ad49eea9a306256c2a3459d90ffc362.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700091_2-2134272/realign && \
touch alignment/FD25700091_2-2134272/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/FD25700091_2-2134272.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/FD25700091_2-2134272.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700091_2-2134272/FD25700091_2-2134272.sorted.realigned.bam
gatk_indel_realigner.FD25700091_2-2134272.9ad49eea9a306256c2a3459d90ffc362.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_102_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_102_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_102_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_103_JOB_ID: gatk_indel_realigner.FD25700092_2-2130593
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700092_2-2130593
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700092_2-2130593.3d52fbd6da5ceddd61d9dad24b0f767c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700092_2-2130593.3d52fbd6da5ceddd61d9dad24b0f767c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700092_2-2130593/realign && \
touch alignment/FD25700092_2-2130593/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/FD25700092_2-2130593.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/FD25700092_2-2130593.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700092_2-2130593/FD25700092_2-2130593.sorted.realigned.bam
gatk_indel_realigner.FD25700092_2-2130593.3d52fbd6da5ceddd61d9dad24b0f767c.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_103_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_103_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_103_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_104_JOB_ID: gatk_indel_realigner.FD25700173_2-2130691
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700173_2-2130691
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700173_2-2130691.029aeadb764136f12b4c1220120033ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700173_2-2130691.029aeadb764136f12b4c1220120033ba.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700173_2-2130691/realign && \
touch alignment/FD25700173_2-2130691/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/FD25700173_2-2130691.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/FD25700173_2-2130691.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700173_2-2130691/FD25700173_2-2130691.sorted.realigned.bam
gatk_indel_realigner.FD25700173_2-2130691.029aeadb764136f12b4c1220120033ba.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_104_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_104_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_104_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_105_JOB_ID: gatk_indel_realigner.FD25700189_2-2130762
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700189_2-2130762
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700189_2-2130762.f58de389ecdff98e95a296ba6e2e9857.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700189_2-2130762.f58de389ecdff98e95a296ba6e2e9857.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700189_2-2130762/realign && \
touch alignment/FD25700189_2-2130762/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/FD25700189_2-2130762.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/FD25700189_2-2130762.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700189_2-2130762/FD25700189_2-2130762.sorted.realigned.bam
gatk_indel_realigner.FD25700189_2-2130762.f58de389ecdff98e95a296ba6e2e9857.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_105_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_105_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_105_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_106_JOB_ID: gatk_indel_realigner.FD25700190_2-2130761
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700190_2-2130761
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700190_2-2130761.8126ca64f7ce14cc990b1cd758c4bdec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700190_2-2130761.8126ca64f7ce14cc990b1cd758c4bdec.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700190_2-2130761/realign && \
touch alignment/FD25700190_2-2130761/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/FD25700190_2-2130761.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/FD25700190_2-2130761.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700190_2-2130761/FD25700190_2-2130761.sorted.realigned.bam
gatk_indel_realigner.FD25700190_2-2130761.8126ca64f7ce14cc990b1cd758c4bdec.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_106_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_106_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_106_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_107_JOB_ID: gatk_indel_realigner.FD25700191_2-2130760
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700191_2-2130760
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700191_2-2130760.11345be9b0fab91273576f955ab7b3dd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700191_2-2130760.11345be9b0fab91273576f955ab7b3dd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700191_2-2130760/realign && \
touch alignment/FD25700191_2-2130760/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/FD25700191_2-2130760.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/FD25700191_2-2130760.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700191_2-2130760/FD25700191_2-2130760.sorted.realigned.bam
gatk_indel_realigner.FD25700191_2-2130760.11345be9b0fab91273576f955ab7b3dd.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_107_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_107_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_107_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_108_JOB_ID: gatk_indel_realigner.FD25700192_2-2130759
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700192_2-2130759
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700192_2-2130759.57ad7489e27ef879b1c192e84e079f05.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700192_2-2130759.57ad7489e27ef879b1c192e84e079f05.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700192_2-2130759/realign && \
touch alignment/FD25700192_2-2130759/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/FD25700192_2-2130759.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/FD25700192_2-2130759.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700192_2-2130759/FD25700192_2-2130759.sorted.realigned.bam
gatk_indel_realigner.FD25700192_2-2130759.57ad7489e27ef879b1c192e84e079f05.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_108_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_108_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_108_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_109_JOB_ID: gatk_indel_realigner.FD25700194_2-2130765
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700194_2-2130765
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700194_2-2130765.8e6ea5c233c28b8560e94f70a2a7031c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700194_2-2130765.8e6ea5c233c28b8560e94f70a2a7031c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700194_2-2130765/realign && \
touch alignment/FD25700194_2-2130765/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/FD25700194_2-2130765.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/FD25700194_2-2130765.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700194_2-2130765/FD25700194_2-2130765.sorted.realigned.bam
gatk_indel_realigner.FD25700194_2-2130765.8e6ea5c233c28b8560e94f70a2a7031c.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_109_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_109_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_109_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_110_JOB_ID: gatk_indel_realigner.FD25700195_2-2130764
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700195_2-2130764
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700195_2-2130764.e451f3928205fd90076879254b7510c7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700195_2-2130764.e451f3928205fd90076879254b7510c7.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700195_2-2130764/realign && \
touch alignment/FD25700195_2-2130764/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/FD25700195_2-2130764.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/FD25700195_2-2130764.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700195_2-2130764/FD25700195_2-2130764.sorted.realigned.bam
gatk_indel_realigner.FD25700195_2-2130764.e451f3928205fd90076879254b7510c7.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_110_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_110_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_110_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_111_JOB_ID: gatk_indel_realigner.FD25700196_2-2130763
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.FD25700196_2-2130763
JOB_DEPENDENCIES=
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.FD25700196_2-2130763.37c1711dc03c4a27549fc5a0943cbac9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.FD25700196_2-2130763.37c1711dc03c4a27549fc5a0943cbac9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/FD25700196_2-2130763/realign && \
touch alignment/FD25700196_2-2130763/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/FD25700196_2-2130763.sorted.bam \
   \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/realign/all.intervals \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx12G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
   \
  --input_file /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/FD25700196_2-2130763.sorted.bam \
   \
  --targetIntervals /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/realign/all.intervals \
   \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/realign/all.bam \
  --maxReadsInMemory 750000 && \
ln -s -f \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/realign/all.bam \
  /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25700196_2-2130763/FD25700196_2-2130763.sorted.realigned.bam
gatk_indel_realigner.FD25700196_2-2130763.37c1711dc03c4a27549fc5a0943cbac9.mugqic.done
chmod 755 $COMMAND
gatk_indel_realigner_111_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem 16G -c 4 -N 1 -q centos7  | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_111_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$gatk_indel_realigner_111_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'10.80.49.1-DnaSeq-FD25671004_2-2130693.FD25671004_2-2130693,FD25671017_2-2134265.FD25671017_2-2134265,FD25671018_2-2134267.FD25671018_2-2134267,FD25671019_2-2134269.FD25671019_2-2134269,FD25671020_2-2134271.FD25671020_2-2134271,FD25671021_2-2134273.FD25671021_2-2134273,FD25671022_2-2134275.FD25671022_2-2134275,FD25671023_2-2134277.FD25671023_2-2134277,FD25671024_2-2134279.FD25671024_2-2134279,FD25671032_2-2130667.FD25671032_2-2130667,FD25672250_2-2130755.FD25672250_2-2130755,FD25672251_2-2130753.FD25672251_2-2130753,FD25672252_2-2130751.FD25672252_2-2130751,FD25672253_2-2130699.FD25672253_2-2130699,FD25672254_2-2130697.FD25672254_2-2130697,FD25672256_2-2130694.FD25672256_2-2130694,FD25672257_2-2130757.FD25672257_2-2130757,FD25672361_2-2130695.FD25672361_2-2130695,FD25672362_2-2130696.FD25672362_2-2130696,FD25672363_2-2130698.FD25672363_2-2130698,FD25672364_2-2130700.FD25672364_2-2130700,FD25672365_2-2130752.FD25672365_2-2130752,FD25672366_2-2130754.FD25672366_2-2130754,FD25672367_2-2130756.FD25672367_2-2130756,FD25672368_2-2130758.FD25672368_2-2130758,FD25673014_2-2130766.FD25673014_2-2130766,FD25699810_2-2130688.FD25699810_2-2130688,FD25699821_2-2130573.FD25699821_2-2130573,FD25699822_2-2130576.FD25699822_2-2130576,FD25699823_2-2130579.FD25699823_2-2130579,FD25699829_2-2130591.FD25699829_2-2130591,FD25699830_2-2130594.FD25699830_2-2130594,FD25699831_2-2130582.FD25699831_2-2130582,FD25699832_2-2130585.FD25699832_2-2130585,FD25699833_2-2130588.FD25699833_2-2130588,FD25699894_2-2130670.FD25699894_2-2130670,FD25699895_2-2130673.FD25699895_2-2130673,FD25699897_2-2130685.FD25699897_2-2130685,FD25699898_2-2130676.FD25699898_2-2130676,FD25699899_2-2130679.FD25699899_2-2130679,FD25699900_2-2130682.FD25699900_2-2130682,FD25699901_2-2134255.FD25699901_2-2134255,FD25699902_2-2134257.FD25699902_2-2134257,FD25699903_2-2134259.FD25699903_2-2134259,FD25699904_2-2134261.FD25699904_2-2134261,FD25699905_2-2134099.FD25699905_2-2134099,FD25699906_2-2134251.FD25699906_2-2134251,FD25699907_2-2134253.FD25699907_2-2134253,FD25699909_2-2130590.FD25699909_2-2130590,FD25699910_2-2130581.FD25699910_2-2130581,FD25699911_2-2130584.FD25699911_2-2130584,FD25699912_2-2130587.FD25699912_2-2130587,FD25699913_2-2130578.FD25699913_2-2130578,FD25699914_2-2130572.FD25699914_2-2130572,FD25699915_2-2130575.FD25699915_2-2130575,FD25699916_2-2130669.FD25699916_2-2130669,FD25699917_2-2130663.FD25699917_2-2130663,FD25699918_2-2130666.FD25699918_2-2130666,FD25699919_2-2130657.FD25699919_2-2130657,FD25699920_2-2130660.FD25699920_2-2130660,FD25699921_2-2130598.FD25699921_2-2130598,FD25699922_2-2130651.FD25699922_2-2130651,FD25699923_2-2130654.FD25699923_2-2130654,FD25699941_2-2130674.FD25699941_2-2130674,FD25699942_2-2130677.FD25699942_2-2130677,FD25699952_2-2130689.FD25699952_2-2130689,FD25699953_2-2130680.FD25699953_2-2130680,FD25699954_2-2130683.FD25699954_2-2130683,FD25699955_2-2130686.FD25699955_2-2130686,FD25699956_2-2130671.FD25699956_2-2130671,FD25699964_2-2130692.FD25699964_2-2130692,FD25699986_2-2130687.FD25699986_2-2130687,FD25699987_2-2130690.FD25699987_2-2130690,FD25699989_2-2130672.FD25699989_2-2130672,FD25699990_2-2130675.FD25699990_2-2130675,FD25699991_2-2130678.FD25699991_2-2130678,FD25699992_2-2130681.FD25699992_2-2130681,FD25699993_2-2130684.FD25699993_2-2130684,FD25700005_2-2130662.FD25700005_2-2130662,FD25700006_2-2130665.FD25700006_2-2130665,FD25700007_2-2130668.FD25700007_2-2130668,FD25700013_2-2134276.FD25700013_2-2134276,FD25700014_2-2134278.FD25700014_2-2134278,FD25700015_2-2134280.FD25700015_2-2134280,FD25700016_2-2130597.FD25700016_2-2130597,FD25700017_2-2130600.FD25700017_2-2130600,FD25700018_2-2130653.FD25700018_2-2130653,FD25700019_2-2130656.FD25700019_2-2130656,FD25700020_2-2130659.FD25700020_2-2130659,FD25700042_2-2130574.FD25700042_2-2130574,FD25700043_2-2130577.FD25700043_2-2130577,FD25700044_2-2130580.FD25700044_2-2130580,FD25700045_2-2130592.FD25700045_2-2130592,FD25700046_2-2130595.FD25700046_2-2130595,FD25700047_2-2130583.FD25700047_2-2130583,FD25700048_2-2130586.FD25700048_2-2130586,FD25700049_2-2130589.FD25700049_2-2130589,FD25700087_2-2134274.FD25700087_2-2134274,FD25700088_2-2134266.FD25700088_2-2134266,FD25700089_2-2134268.FD25700089_2-2134268,FD25700090_2-2134270.FD25700090_2-2134270,FD25700091_2-2134272.FD25700091_2-2134272,FD25700092_2-2130593.FD25700092_2-2130593,FD25700173_2-2130691.FD25700173_2-2130691,FD25700189_2-2130762.FD25700189_2-2130762,FD25700190_2-2130761.FD25700190_2-2130761,FD25700191_2-2130760.FD25700191_2-2130760,FD25700192_2-2130759.FD25700192_2-2130759,FD25700194_2-2130765.FD25700194_2-2130765,FD25700195_2-2130764.FD25700195_2-2130764,FD25700196_2-2130763.FD25700196_2-2130763' | md5sum | awk '{ print $1 }')
if test -t 1; then ncolors=$(tput colors); if test -n "$ncolors" && test $ncolors -ge 8; then bold="$(tput bold)"; normal="$(tput sgr0)"; yellow="$(tput setaf 3)"; fi; fi
wget --quiet 'http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=narval1.narval.calcul.quebec&ip=10.80.49.1&pipeline=DnaSeq&steps=bwa_mem_sambamba_sort_sam,sambamba_merge_sam_files,gatk_indel_realigner,sambamba_merge_realigned&samples=111&md5=$LOG_MD5' -O /dev/null || echo "${bold}${yellow}Warning:${normal}${yellow} Genpipes ran successfully but was not send telemetry to mugqic.hpc.mcgill.ca. This error will not affect genpipes jobs you have submitted.${normal}"
