library(methylKit)
library(tidyverse)
library(genomation)

##module load mugqic/R_Bioconductor/4.1.0_3.13

##salloc --time=2:00:0 --account=rrg-bourqueg-ad --mem=64000M

#library(edgeR)
library(methylKit)
library(tidyverse)
library(genomation)

bismark_files <- as.list(list.files(path="Bismark_files",full.names=T))
sample_id <- gsub("\\.readset_sorted.dedup.bismark.cov.gz","",basename(list.files(path="Bismark_files",full.names=T)))
out <- "Differential_Methylation_1st_iteration"
dir.create(out,showWarnings=F)
#meth=methRead(bismark_files,
#           sample.id=as.list(sample_id),
#           assembly="Okis",
#           context="CpG",
#           treatment=c(rep(0,6),rep(1,6)),
#           mincov = 8,
#          pipeline = "bismarkCoverage",
#          dbtype = "tabix"
#           )
#
#save(meth,file=file.path(out,"meth_object.rda"))
load(file=file.path(out,"meth_object.rda"))
### testing with smalll dataset
#
##meth=methRead(bismark_files[1:4],
##           sample.id=as.list(sample_id)[1:4],
##           assembly="Okis",
##           context="CpG",
##           treatment=c(0,0,1,1),
##           mincov = 8,
##          pipeline = "bismarkCoverage",
##          dbtype = "tabix"
##           )
#


##### QC ######################
#names(meth) <- sample_id
#map2(meth,names(meth),function(x,y) {
#	png(file=file.path(out,paste0(y,"_Coverage.png")))
#		getCoverageStats(x,plot=T)
#	dev.off()	
#})
#
#
#
#map2(meth,names(meth),function(x,y) {
#        png(file=file.path(out,paste0(y,"_Methylation_distribution.png")))
#                getMethylationStats(x,plot=T)
#        dev.off()       
#})
#
#
####### Preprocessing ####################
#filtered.meth=filterByCoverage(meth,lo.count=10,lo.perc=NULL,
#                                hi.count=NULL,hi.perc=99.9)
#
#
#save(filtered.meth,file=file.path(out,"filtered_meth_obj.rda"))
#normed.meth <- normalizeCoverage(filtered.meth)
#
#save(normed.meth,file=file.path(out,"normed_meth_obj.rda"))
#
#load(file=file.path(out,"normed_meth_obj.rda"))
#unite.meth=methylKit::unite(normed.meth)
#
#save(unite.meth,file=file.path(out,"unite_meth_object.rda"))
#unite.meth
#
#png(file=file.path(out,"Correlation_plot.png"))
#	getCorrelation(unite.meth,plot=TRUE)
#dev.off()
#
#
#png(file=file.path(out,"SampleClustering.png"))
#	clusterSamples(unite.meth, dist="correlation", method="ward", plot=TRUE)
#dev.off()
#
#
#png(file=file.path(out,"PCA_plot_all_regions.png"))
#	PCASamples(unite.meth)
#dev.off()


load(file=file.path(out,"unite_meth_object.rda"))

### added later
# better pca plot
meta <- read.csv("Metadata.csv")
pca <- PCASamples(unite.meth,obj.return=T)

pca_plot_df <- as.data.frame(pca$x) %>% 
rownames_to_column("Sample") %>%
dplyr::select(Sample:PC2) %>%
inner_join(.,meta,by="Sample") %>%
dplyr::select(-File) %>% 
mutate(timepoint = as.factor(timepoint))


var <- round(100*(pca$sdev)^2/sum(pca$sdev^2),2)
library(ggrepel)

p <- ggplot(pca_plot_df,aes(PC1,PC2)) +
geom_point(aes(color=timepoint,shape=genotype),size=3) +
labs(x=paste0("PC1: ",var[1],"% variance"),y=paste0("PC2: ",var[2], "% variance")) +
scale_shape_manual(values=c(16,2)) +
geom_text_repel(data=filter(pca_plot_df,PC1>(1000)),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_all.png"))


######################### Preparing annotation #####################
##library(GenomicFeatures)
##ibrary("AnnotationDbi")
##
##
##
##gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
##                                           package = "methylKit"))
##
##
##
##gtf <- "../../genomes/Okis_V2_genomic.gtf"
##gtf_gr <- makeTxDbFromGFF(gtf)
##Import genomic features from the file as a GRanges object ... OK
##Prepare the 'metadata' data frame ... OK
##Make the TxDb object ... OK
##Warning messages:
##1: In .get_cds_IDX(mcols0$type, mcols0$phase) :
##  The "phase" metadata column contains non-NA values for features of type
##  stop_codon. This information was ignored.
##2: In makeTxDbFromGRanges(gr, metadata = metadata) :
##  The following transcripts were dropped because their exon ranks could
##  not be inferred (either because the exons are not on the same
##  chromosome/strand or because they are not separated by introns):
##  unknown_transcript_1
##3: In makeTxDbFromGRanges(gr, metadata = metadata) :
##  The following transcripts were dropped because no genomic ranges could
##  be found for them and their ranges could not be inferred from their
##  exons either (because they have them on both strands):
##  unknown_transcript_1
##
##save(gtf_gr,file=file.path(out,"txdb.rda"))
###exons <- GenomicFeatures::exonsBy(gtf_gr,use.names=T)
####exons <- mapply(function(x,n) {x$name <- n; return(x)},exons,names(exons))
##
##
##exons_gr <- unique(exons(gtf_gr,columns=c("exon_rank","exon_id","gene_id")))
### promoters - upstream 2000bp downstream 200bp
##promoters_gr <- unique(promoters(gtf_gr,columns = c("tx_id","tx_name","gene_id")))
##introns_gr <- intronsByTranscript(gtf_gr,use.names=T)
##
##introns_gr <- unlist(introns_gr,use.names=T)
##introns_gr <- unique(introns_gr)
##introns_gr$txid <- names(introns_gr)
##names(introns_gr) <- 1:length(introns_gr)
##genes_gr <- unique(genes(gtf_gr))
##
##save(introns_gr,genes_gr,exons_gr,promoters_gr,file="annotation_granges.rda")
##
##
##load(file="annotation_granges.rda")
##
##
###### DM analysis
##files <- list.files(path="Bismark_files",full.names=T)
##sample_id <- gsub("\\.readset_sorted.dedup.bismark.cov.gz","",basename(files))
group <- gsub("DNA.","",sample_id)
genotype <- c(rep("NT",6),rep("TF",6))
timepoint <- rep(rep(c(0,2,3),each=2),2)
targets <- data.frame(Sample=sample_id,
                      group=group,
                      genotype=genotype,
                      timepoint=timepoint,
                      File=unlist(bismark_files))
#write.csv(targets,file="Metadata.csv",row.names=F)

comparisons <- list(
TF_Day0_vs_TF_Starve28 = c("TF2","TF0"),
 TF_Day0_vs_TF_Refed41 = c("TF3", "TF0"),
 TF_Starve28_vs_TF_Refed41 = c("TF3","TF2"),

 NT_Day0_vs_NT_Starve28 = c("NT2","NT0"),
 NT_Day0_vs_NT_Refed41 = c("NT3","NT0"),
 NT_Starve28_vs_NT_Refed41 = c("NT3","NT2"),

 NT_Day0_vs_TF_Day0 = c("TF0","NT0"),
 NT_Starve28_vs_TF_Starve28 = c("TF2", "NT2"),
NT_Refed41_vs_TF_Refed41 =c("TF3","NT3")
)



#### for overall
#meth_objs_for_dm <- map(comparisons,function(c) {
#treat <- targets[targets$group == c[1],"Sample"]
#control <- targets[targets$group == c[2],"Sample"]
#
#reorganize(unite.meth,sample.ids=c(treat,control),treatment = c(1,1,0,0))
#
#})
#
#save(meth_objs_for_dm ,file=file.path(out,"objects_for_dm.rda"))
#
#
#Diffs <- map2(meth_objs_for_dm,names(meth_objs_for_dm),function(x,n) {
#	diff_stats <- calculateDiffMeth(x)
#	save(diff_stats,file=file.path(out,paste0(n,"_DM_stat_objects.rda")))
#
#	diff20p <- getMethylDiff(diff_stats,difference=20,qvalue=0.05)
#	save(diff20p,file=file.path(out,paste0(n,"_DM_20p.rda")))
#	
#	 diff20p_gr <- as(diff20p,"GRanges")
#	write.csv(as.data.frame(diff20p_gr),file=file.path(out,paste0(n,"_global_DMtable.csv")),row.names=F)
#
#	return(diff20p)
#})
#
#
#
###### By Genomic regions
load(file="annotation_granges.rda")

#map2(list(promoters_gr,genes_gr,exons_gr,introns_gr),c("promoters_gr","genes_gr","exons_gr","introns_gr"),function(x,n) {
map2(list(introns_gr),list("introns_gr"),function(x,n) {
meth <- regionCounts(meth,x,save.db=F)
	filtered.meth=filterByCoverage(meth,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9,suffix = n)


	save(filtered.meth,file=file.path(out,paste0(n,"_filtered_meth_obj.rda")))
	normed.meth <- normalizeCoverage(filtered.meth,suffix=n)
:
	save(normed.meth,file=file.path(out,paste0(n,"_normed_meth_obj.rda")))
	unite.meth=methylKit::unite(normed.meth,suffix = n)

	save(unite.meth,file=file.path(out,paste0(n,"_unite_meth_object.rda")))
	unite.meth

	png(file=file.path(out,paste0(n,"_Correlation_plot.png")))
		getCorrelation(unite.meth,plot=TRUE)
	dev.off()


	png(file=file.path(out,paste0(n,"_SampleClustering.png")))
        	clusterSamples(unite.meth, dist="correlation", method="ward", plot=TRUE)
	dev.off()


	png(file=file.path(out,paste0(n,"_PCA_plots.png")))
        	PCASamples(unite.meth)
	dev.off()


	meth_objs_for_dm <- map(comparisons,function(c) {
	treat <- targets[targets$group == c[1],"Sample"]
	control <- targets[targets$group == c[2],"Sample"]

	reorganize(unite.meth,sample.ids=c(treat,control),treatment = c(1,1,0,0))

	})

	save(meth_objs_for_dm ,file=file.path(out,paste0(n,"_objects_for_dm.rda")))


	Diffs <- map2(meth_objs_for_dm,names(meth_objs_for_dm),function(x,c) {
        	diff_stats <- calculateDiffMeth(x)
        	save(diff_stats,file=file.path(out,paste(c,n,"DM_stat_objects.rda",sep="_")))

        	diff20p <- getMethylDiff(diff_stats,difference=20,qvalue=0.05)
        	save(diff20p,file=file.path(out,paste(c,n,"DM_20p.rda",sep="_")))

		diff20p_gr <- as(diff20p,"GRanges")
        	write.csv(as.data.frame(diff20p_gr),file=file.path(out,paste(c,n,"DMtable.csv",sep="_")),row.names=F)

       		 return(diff20p)
	})
})


### redoing the PCA plots
load(file=file.path(out,"promoters_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,obj.return=T)

var <- round(100*(pca$sdev)^2/sum(pca$sdev^2),2)

## promoters
pca_plot_df <- as.data.frame(pca$x) %>%
rownames_to_column("Sample") %>%
dplyr::select(Sample:PC2) %>%
inner_join(.,meta,by="Sample") %>%
dplyr::select(-File) %>%
mutate(timepoint = as.factor(timepoint))


p <- ggplot(pca_plot_df,aes(PC1,PC2)) +
geom_point(aes(color=timepoint,shape=genotype),size=3) +
labs(x=paste0("PC1: ",var[1],"% variance"),y=paste0("PC2: ",var[2],"% variance")) +
scale_shape_manual(values=c(16,2)) +
geom_text_repel(data=filter(pca_plot_df,PC1>200),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_promoters_gr.png"))

### exons
load(file=file.path(out,"exons_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,obj.return=T)


pca_plot_df <- as.data.frame(pca$x) %>%
rownames_to_column("Sample") %>%
dplyr::select(Sample:PC2) %>%
inner_join(.,meta,by="Sample") %>%
dplyr::select(-File) %>%
mutate(timepoint = as.factor(timepoint))


p <- ggplot(pca_plot_df,aes(PC1,PC2)) +
geom_point(aes(color=timepoint,shape=genotype),size=3) +
labs(x="PC1",y="PC2") +
scale_shape_manual(values=c(16,2)) +
geom_text_repel(data=filter(pca_plot_df,PC1>200),aes(label=Sample)) +
theme_bw()


ggsave(p,file="PCA_exons.png")


## introns
load(file=file.path(out,"introns_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,obj.return=T)


pca_plot_df <- as.data.frame(pca$x) %>%
rownames_to_column("Sample") %>%
dplyr::select(Sample:PC2) %>%
inner_join(.,meta,by="Sample") %>%
dplyr::select(-File) %>%
mutate(timepoint = as.factor(timepoint))


p <- ggplot(pca_plot_df,aes(PC1,PC2)) +
geom_point(aes(color=timepoint,shape=genotype),size=3) +
labs(x="PC1",y="PC2") +
scale_shape_manual(values=c(16,2)) +
geom_text_repel(data=filter(pca_plot_df,PC1>400),aes(label=Sample)) +
theme_bw()


ggsave(p,file="PCA_introns.png")

## genebodies
load(file=file.path(out,"genes_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,obj.return=T)


pca_plot_df <- as.data.frame(pca$x) %>%
rownames_to_column("Sample") %>%
dplyr::select(Sample:PC2) %>%
inner_join(.,meta,by="Sample") %>%
dplyr::select(-File) %>%
mutate(timepoint = as.factor(timepoint))


p <- ggplot(pca_plot_df,aes(PC1,PC2)) +
geom_point(aes(color=timepoint,shape=genotype),size=3) +
labs(x="PC1",y="PC2") +
scale_shape_manual(values=c(16,2)) +
geom_text_repel(data=filter(pca_plot_df,PC1>200),aes(label=Sample)) +
theme_bw()


ggsave(p,file="PCA_genebodies.png")









##targets
##
##yall <- readBismark2DGE(targets$File, sample.names=targets$Sample)
##ytest <- readBismark2DGE(targets$File[1:4], sample.names=targets$Sample[1:4])
#Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
#
#Me <- yall$counts[, Methylation=="Me"]
#Un <- yall$counts[, Methylation=="Un"]
#Coverage <- Me + Un
#HasCoverage <- rowSums(Coverage >= 8) == 2
#
#table(HasCoverage)
##HasCoverage
##   FALSE     TRUE
##43684593 27175110
#
#
#
#HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0
#
#table(HasCoverage,HasBoth)
#
#
#y <- yall[HasCoverage & HasBoth,,keep.lib.sizes=F]
#
#
#TotalLibSize <- y$samples$lib.size[Methylation=="Me"] + y$samples$lib.size[Methylation=="Un"]
#y$samples$lib.size <- rep(TotalLibSize, each=2)
#
#
#
#Me <- y$counts[, Methylation=="Me"]
#Un <- y$counts[, Methylation=="Un"]
#M <- log2(Me + 2) - log2(Un + 2)
#colnames(M) <- targets$Sample
##colnames(M) <- targets$Sample[1:4]
#
#
#
#mds <- plotMDS(M)
#
### create better plot
#mds_df <- data.frame(Sample=colnames(M),
#		    x=mds$x,
#		    y=mds$y,
#		    GT=targets$genotype,
#		    Time=as.factor(targets$timepoint))
#
#
##mds_df <- data.frame(Sample=colnames(M),
##                    x=mds$x,
##                    y=mds$y,
##                    GT=targets$genotype[1:4],
##                    Time=as.factor(targets$timepoint[1:4]))
##
#
#ggplot(mds_df,aes(x,y))+
#geom_point(aes(color=GT,shape=Time)) +
#xlab("Leading LogFC dim 1") +
#ylab("Leading LogFC dim 2") +
#scale_color_hue("Genotype") +
#scale_shape(breaks=c(0,2,3),label=c("Day 0","Starve 28","Refed 41")) +
#theme_bw()
#
#ggsave(file = file.path(out,"MDS.pdf"))
#
#
##pdf(file="MDS.pdf")
##plotMDS(M, col=rep(1:3, 2), main="M-values")
##dev.off()
#
#
##designSL <- model.matrix(~0+group, data=targets[1:4,])
#designSL <- model.matrix(~0+group, data=targets)
#design <- modelMatrixMeth(designSL)
#y <- estimateDisp(y, design=design, trend="none")
#fit <- glmFit(y, design)
#
#save(y,file=file.path(out,"edgeR_object.rda"))
##y1 <- y[y$genes$Chr == "NC_034184.2",]
#
### make contrasts
#contr <- makeContrasts(
# TF_Day0_vs_TF_Starve28 = groupTF2 - groupTF0,
# TF_Day0_vs_TF_Refed41 = groupTF3 - groupTF0,
# TF_Starve28_vs_TF_Refed41 = groupTF3 - groupTF2,
#
# NT_Day0_vs_NT_Starve28 = groupNT2 - groupNT0,
# NT_Day0_vs_NT_Refed41 = groupNT3 - groupNT0,
# NT_Starve28_vs_NT_Refed41 = groupNT3 - groupNT2,
#
# NT_Day0_vs_TF_Day0 = groupTF0 - groupNT0,
# NT_Starve28_vs_TF_Starve28 = groupTF2 - groupNT2,
# NT_Refed41_vs_TF_Refed41 = groupTF3 - groupNT3,
# levels=design)
#
#lrts <- apply(contr,2,function(x) {glmLRT(fit, contrast=x})
#save(lrts,file=file.path(out,"lrts.rda"))
#
#names(lrts) <- colnames(contr)
#
#file_names <- file.path(out,paste0(names(lrt),"_DM_table.csv"))
#
#map2(lrts,filenames,function(x,y) {write.csv(topTags(x,n=Inf),file=y,row.names=F)})
##contr <- makeContrasts(NF_Day0_vs_NF_Starve28 = groupNT2 - groupNT0,levels=design)
