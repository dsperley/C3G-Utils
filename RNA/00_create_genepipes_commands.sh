#!/bin/bash


rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.cedar.ini ./Susta.ini -r ./ReadSet.txt -d  ./design.txt -s 1-17 --no-json -t stringtie -o ./genpipes -g 01_rnaseq_commands_s1_17.sh

#rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.cedar.ini ./Delphinapterus_leucas.ASM228892v3.ini --report -r ./ReadSet.txt -d  ./design.txt -s 1-17 --no-json -t stringtie -o ./genpipes > 01_rnaseq_commands_report.sh
