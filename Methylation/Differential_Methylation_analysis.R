
##module load mugqic/R_Bioconductor/4.1.0_3.13

##salloc --time=2:00:0 --account=rrg-bourqueg-ad --mem=64000M

#library(edgeR)
library(methylKit)
library(tidyverse)
library(genomation)



args <- commandArgs(trailingOnly=T)
## script inputDir outDir
dir <- args[1]
out <- args[2]
path <- "/home/dsperley/projects/rrg-bourqueg-ad/C3G/projects/BDelvin_Okis_C3G-RD-542/WGBS/genpipes_output/methylkit/inputs"
bismark_files <- as.list(list.files(path=path,pattern="input",full.names=T))
sample_id <- gsub("\\.readset_sorted.dedup.map.input","",basename(list.files(path=path,pattern = "input")))

dir.create(out,showWarnings=F)
meta <- read.csv("Metadata.csv")

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

#save(meth_objs_for_dm ,file=file.path(out,"objects_for_dm.rda"))
 load(file.path(dir,"objects_for_dm.rda"))

Diffs <- map2(meth_objs_for_dm,names(meth_objs_for_dm),function(x,n) {
	diff_stats <- calculateDiffMeth(x,overdispersion="MN",mc.cores=8)
	save(diff_stats,file=file.path(out,paste0(n,"_DM_stat_objects.rda")))

	diff20p <- getMethylDiff(diff_stats,difference=20,qvalue=0.05)
	save(diff20p,file=file.path(out,paste0(n,"_DM_20p.rda")))
	
	 diff20p_gr <- as(diff20p,"GRanges")
	write.csv(as.data.frame(diff20p_gr),file=file.path(out,paste0(n,"_global_DMtable.csv")),row.names=F)

	return(diff20p)
})



###### By Genomic regions
load(file="annotation_granges.rda")

map2(list(promoters_gr,genes_gr,exons_gr,introns_gr),c("promoters_gr","genes_gr","exons_gr","introns_gr"),function(x,n) {
#map2(list(introns_gr),list("introns_gr"),function(x,n) {
#meth <- regionCounts(meth,x,save.db=F)

	load(file=file.path(dir,paste0(n,"_objects_for_dm.rda")))

	Diffs <- map2(meth_objs_for_dm,names(meth_objs_for_dm),function(x,c) {
        	diff_stats <- calculateDiffMeth(x,overdispersion="MN",mc.cores=8)
        	save(diff_stats,file=file.path(out,paste(c,n,"DM_stat_objects.rda",sep="_")))

        	diff20p <- getMethylDiff(diff_stats,difference=20,qvalue=0.05)
        	save(diff20p,file=file.path(out,paste(c,n,"DM_20p.rda",sep="_")))

		diff20p_gr <- as(diff20p,"GRanges")
        	write.csv(as.data.frame(diff20p_gr),file=file.path(out,paste(c,n,"DMtable.csv",sep="_")),row.names=F)

       		 return(diff20p)
	})
})


