#!/usr/bin/env Rscript
library(stringr)
library(purrr) 



args <- commandArgs(trailingOnly=TRUE)

arf <- args[1]
out <- args[2]
arf <- read.table(arf,header=F)

arf <- arf[!duplicated(arf$V1),]
arf$sample <- str_split(arf$V1,"_",simplify=T)[,1]
arf$n_reads <- as.numeric(str_match(arf$V1,"_x(\\d+)")[,2]) 

arf_by_sample <- split(arf,arf$sample)  
total_per_sample <- map(arf_by_sample,~sum(.x$n_reads)) 

virus_per_sample <- map(arf_by_sample,function(x) sum(x[x$V6=="KT819632.1","n_reads"])) 
stats <- map2_dfr(total_per_sample,virus_per_sample,function(x,y)  {total=x;virus=y; perc.virus=round((y*100/x),2);perc.host=100-perc.virus; data.frame(total_aligned=total,virus_reads=virus,percentage_virus=perc.virus,percentage_host=perc.host)},.id="Sample")

write.table(stats,file=out,sep="\t",row.names=F,quote=F)
