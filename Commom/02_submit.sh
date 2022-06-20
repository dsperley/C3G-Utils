#!/bin/bash

#$MUGQIC_PIPELINES_HOME/utils/chunk_genpipes.sh 01_rnaseq_commands_s1_17.sh 01_rnaseq_commands_s1_17_chunks

$MUGQIC_PIPELINES_HOME/utils/submit_genpipes 01_rnaseq_commands_s1_17_chunks -n 20

