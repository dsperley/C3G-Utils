#!/bin/bash

M_FOLDER=dnaseq_s8_to_end

$MUGQIC_PIPELINES_HOME/utils/chunk_genpipes.sh 05_dnaseq_s8_to_end.sh $M_FOLDER

$MUGQIC_PIPELINES_HOME/utils/submit_genpipes  $M_FOLDER
