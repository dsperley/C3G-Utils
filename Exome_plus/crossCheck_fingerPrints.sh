sbatch --mail-user=eveleigh.rjm@gmail.com --mail-type=end,fail --account=${RAP_ID} --job-name="fingerprint.${ID}" -n 1 -c 10 --mem=40G --time=71:00:00 \
 --wrap="module purge && module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/4.2.2.0 && \
 gatk \
 --java-options \"-Djava.io.tmpdir=\${TMPDIR:=/tmp} -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx36G\" \
  CrosscheckFingerprints --NUM_THREADS 6 --EXIT_CODE_WHEN_MISMATCH 0 \
  --VALIDATION_STRINGENCY SILENT \
  --TMP_DIR \${TMPDIR:=/tmp} \
  --REFERENCE_SEQUENCE \$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --HAPLOTYPE_MAP \$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.fingerprint.map \
  --LOD_THRESHOLD 3.0 \
  --OUTPUT PAD.sample.fingerprint \
  --MATRIX_OUTPUT PAD.sample.fingerprint.matrix \
  --INPUT bam_list.tsv"