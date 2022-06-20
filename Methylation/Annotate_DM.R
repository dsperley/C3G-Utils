## module load mugqic/R_Bioconductor/4.1.0_3.13


args <- commandArgs(trailingOnly=T)
dir <- args[1]
out <- args[1]
library(GenomicRanges)
library("AnnotationDbi")
library(tidyverse)
library(ChIPseeker)
load(file="annotation_granges.rda")
gtf_gr <- loadDb(file = "Differential_Methylation_1st_iteration/txdb2.rda")
dir <- out <- "Differential_Methylation_1.5_iteration_overdispersion"
dir.create(out,showWarnings=F)

dm_files <- list.files(path=dir,pattern= "DMtable.csv",full.names=T)
dm_files
dm_res <- lapply(dm_files,read.csv)

names(dm_res) <- gsub("_DMtable.csv","",basename(dm_files))

dm_gr <- lapply(dm_res,makeGRangesFromDataFrame,keep.extra.columns=T)

#as.data.frame(mcols(test_add))

annot <- function(x,y) {
hits <- findOverlaps(x,y,type="equal")
mcols(x) <- c(mcols(x),mcols(y[subjectHits(hits),]))

dat <- as.data.frame(x)

return(dat)
}


## annotating regions
DM_by_region <- list(promoters=dm_gr[grep("promoter",names(dm_gr))],
		     exons=dm_gr[grep("exons",names(dm_gr))],
	    	     genes=dm_gr[grep("genes",names(dm_gr))],
		     introns=dm_gr[grep("introns",names(dm_gr))]
		    )


## add gene names to intron gr

tx_gene <- AnnotationDbi::select(gtf_gr,introns_gr$txid,columns=c("TXNAME","GENEID"),keytype="TXNAME")
tx_gene <- tx_gene[match(introns_gr$txid,tx_gene$TXNAME),]
introns_gr$GENEID <- tx_gene$GENEID
### fix formatting

dm_annotated <- mapply(function(x,y){lapply(x,annot,y)},x=DM_by_region,y=list(promoters_gr,exons_gr,genes_gr,introns_gr))
names(dm_annotated) <-  unlist(lapply(DM_by_region,names))

## fix promoters and exons
dm_annotated[grep("exons",names(dm_annotated))] <- map(dm_annotated[grep("exons",names(dm_annotated))],function(x) {
x$exon_rank <- sapply(x$exon_rank,function(i) paste(unique(i),collapse=";"))
x$gene_id <- sapply(x$gene_id,function(i) paste(unique(i),collapse=";"))
return(x)
})


dm_annotated[grep("promoters",names(dm_annotated))] <- map(dm_annotated[grep("promoters",names(dm_annotated))],function(x) {
x$gene_id <- unlist(x$gene_id)
return(x)
})


file_names <- file.path(out,paste0(names(dm_annotated),"_DMtable_annotated.txt"))
mapply(function(x,y) {write.table(x,file=y,row.names=F,sep="\t")},x=dm_annotated,y=file_names)


#### annotate global results
global_gr <- dm_gr[grep("global",names(dm_gr))]
names(global_gr) <- gsub("_global","",names(global_gr))

#test <- global_gr[[1]]

#annot <- annotatePeak(test,tssRegion=c(-2000,200),TxDb=gtf_gr)

annotated <- map2(global_gr,names(global_gr),function(x,n) {
	a <- annotatePeak(x,tssRegion=c(-2000,200),TxDb=gtf_gr)
	
	png(file=file.path(out,paste0(n,"_DM_annotation.png")))
		plotAnnoPie(a)
	dev.off()

	df <- as.data.frame(a,stringsAsFactors=F)
	df2 <- dplyr::select(df,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand)
	return(df2)
	})


map2(annotated,names(annotated),function(x,n) {write.table(x,file=file.path(out,paste0(n,"_global_DMtable_annotated.txt")),row.names=F,sep="\t",quote=F)})	
