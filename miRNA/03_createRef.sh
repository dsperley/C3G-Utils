
module load mugqic/python/3.7.3 mugqic/bowtie/1.2.2 

PROJECT="/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/Pearson_canine_miRNA_R001401"

#cat /cvmfs/soft.mugqic/CentOS6/genomes/species/Canis_familiaris.CanFam3.1/genome/Canis_familiaris.CanFam3.1.fa ${PROJECT}/Ref/KT819632_1.fa > $PROJECT/Ref/Canis_familiaris.CanFam3.1_Canid_herpesvirus.fa

bowtie-build $PROJECT/Ref/Canis_familiaris.CanFam3.1_Canid_herpesvirus.fa $PROJECT/Ref/CamFam3.1_Canid_virus
