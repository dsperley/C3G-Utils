#!/bin/bash

 awk 'BEGIN{FS = "\t";OFS= "\t"} {print $38,$21}' CHV-1_MDCK_microRNAs_NovaSeqReadSet_2021-09-21.txt | tail -n +2 > Adapters.txt 

sed 's/ dna.*//'g Canis_familiaris.CanFam3.1_Canid_herpesvirus.fa | sed 's/ Canid.*//g' > Canis_familiaris.CanFam3.1_Canid_herpesvirus_no_white_spaces.fa  

sed 's/ /_/g' CanFam3.1_hairpins.fa > CanFam3.1_hairpins_no_white_space.fa 
sed 's/ /_/g' CanFam3.1_mature.fa > CanFam3.1_mature_no_white_space.fa 
