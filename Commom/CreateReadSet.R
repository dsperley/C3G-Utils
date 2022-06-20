#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
## CreateReadSet.r in path out


files <- read.csv(args[1])

path <- args[2]
readset <- data.frame(Sample=files$Name,
                      Readset=paste(files$Name,files$Run,files$Region,sep="."),
                      Library=files$Library.Name,
                      RunType=files$Run.Type,
                      Lane=files$Region,
                      Adapter1=files$Adaptor.Read.1..NOTE..Usage.is.bound.by.Illumina.Disclaimer.found.on.Nanuq.Project.Page.,
                      Adapter2=files$Adaptor.Read.2..NOTE..Usage.is.bound.by.Illumina.Disclaimer.found.on.Nanuq.Project.Page.,
                      QualityOffset=files$Quality.Offset,
                      Bed="",
                      FASTQ1=file.path(path,paste(files$Filename.Prefix,"R1.fastq.gz",sep="_")),
                      FASTQ2=file.path(path,paste(files$Filename.Prefix,"R2.fastq.gz",sep="_")),
                      Bam="")



write.table(readset,file=args[3],row.names=F,sep="\t",quote=F)

