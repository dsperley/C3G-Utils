#!/bin/bash

cat <(echo "Sample;Readset;Library;RunType;Run;Lane;Adapter1;Adapter2;QualityOffset;BED;FASTQ1;FASTQ2;BAM" | tr ';' '\t') \
<(for i in `ls raw_reads/*/*R1*.fastq.gz | grep "Cut_0"` ; \
do R2=`echo $i | sed -e 's#R1#R2#g'` ; \
SID=`echo $i | cut -d/ -f2 | sed -e 's#Sample_##g'` ; \
LID=`echo $i | cut -d/ -f3 | sed -e 's#.*_L##g' | sed -e 's#_R[0-9]*.*##'` ; \
RUN=`echo $i | cut -d/ -f3 | sed -e 's#.*_S##g' |  sed -e 's#_L.*##g'` ; \
LANE=`echo $i | cut -d/ -f3 | cut -d_ -f3` ; \
echo "$SID;$SID;$SID;PAIRED_END;S$RUN;L$LID;AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA;AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG;33;/lustre06/project/6061810/dsperley/E100043108/SureSelectHumanAllExonV7.Target.b38.bed;$i;$R2;" ; \
done | tr ';' '\t'  | sort -k1,1V)

