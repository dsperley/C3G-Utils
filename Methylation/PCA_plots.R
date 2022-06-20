##module load mugqic/R_Bioconductor/4.1.0_3.13

##salloc --time=2:00:0 --account=rrg-bourqueg-ad --mem=64000M

#library(edgeR)
library(methylKit)
library(tidyverse)
library(genomation)

#bismark_files <- as.list(list.files(path="Bismark_files",full.names=T))
#sample_id <- gsub("\\.readset_sorted.dedup.bismark.cov.gz","",basename(list.files(path="Bismark_files",full.names=T)))
dir <- "Differential_Methylation_1.5_iteration"
out <- "PCA_plots"
dir.create(out,showWarnings=F)


load(file=file.path(dir,"unite_meth_object.rda"))


# better pca plot
meta <- read.csv("Metadata.csv")
pca <- PCASamples(unite.meth,sd.threshold=0.95,obj.return=T)

pca_plot_df <- as.data.frame(pca$x) %>% 
rownames_to_column("Sample") %>%
dplyr::select(Sample:PC2) %>%
inner_join(.,meta,by="Sample") %>%
dplyr::select(-File) %>% 
mutate(timepoint = as.factor(timepoint))


var <- round(100*(pca$sdev)^2/sum(pca$sdev^2),2)


p <- ggplot(pca_plot_df,aes(PC1,PC2)) +
geom_point(aes(color=timepoint,shape=genotype),size=3) +
labs(x=paste0("PC1: ",var[1],"% variance"),y=paste0("PC2: ",var[2], "% variance")) +
scale_shape_manual(values=c(16,2)) +
#geom_text_repel(data=filter(pca_plot_df,PC1>(1000)),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_all.png"))



### better PCA plots of genomic regions

load(file=file.path(dir,"promoters_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,sd.threshold=0.95,obj.return=T)

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
#geom_text_repel(data=filter(pca_plot_df,PC1>200),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_promoters_gr.png"))

### exons
load(file=file.path(dir,"exons_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,sd.threshold=0.95,obj.return=T)


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
#geom_text_repel(data=filter(pca_plot_df,PC1>200),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_exons.png"))


## introns
load(file=file.path(dir,"introns_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,sd.threshold=0.95,obj.return=T)


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
#geom_text_repel(data=filter(pca_plot_df,PC1>400),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_introns.png"))

## genebodies
load(file=file.path(dir,"genes_gr_unite_meth_object.rda"))
pca <- PCASamples(unite.meth,sd.threshold=0.95,obj.return=T)


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
#geom_text_repel(data=filter(pca_plot_df,PC1>200),aes(label=Sample)) +
theme_bw()


ggsave(p,file=file.path(out,"PCA_genebodies.png"))









