#!/bin/bash

module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/4.1.8.1 mugqic/R_Bioconductor/3.5.1_3.7 && \
gatk --java-options "-Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx6G" \
ValidateSamFile \
      -I /lustre06/project/6061810/dsperley/Run1/genpipes/alignment/FD25699894_2-2130670/FD25699894_2-2130670.sorted.dup.recal.bam \
      --MODE SUMMARY
