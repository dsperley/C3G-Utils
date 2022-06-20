#!/bin/bash

PROJECT="/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401/"

module unload mugqic/python/2.7.14
module load  python/3.9.6 mugqic/fastqc/0.11.5 mugqic/MultiQC/1.10.1

OUT="${PROJECT}/01_FASTQC"

if [ ! -d "$OUT" ]; then mkdir "$OUT"; fi

#fastqc -o $OUT ${PROJECT}/fastq/NS*gz

multiqc ${OUT}/*zip
