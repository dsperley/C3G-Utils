#!/bin/bash

export RAP_ID=rrg-vmooser
module purge && module load mugqic/python/3.10.4 mugqic/genpipes/4.1.0 && dnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.base.ini $MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.exome.ini $MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.beluga.ini Exome_plus.ini -r readSet.txt -s 3-6 -o genpipes -g 03_dnaseq_s3-6b.sh
