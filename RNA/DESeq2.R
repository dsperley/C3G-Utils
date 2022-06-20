

################ Setup ####################################
library("DESeq2")
library("gplots")
library("RColorBrewer") 
library("ggplot2")
library("tidyverse")
library("pheatmap") 
library("biomaRt")
library("ashr")
library("readxl")
library("gprofiler2")
library("goseq")


mytheme <-  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
theme_set(mytheme)

#setwd("/Volumes/GoogleDrive/My Drive/Projects/Tanya_Beluga_RNAseq_R008586/")
source("./RNA_DE_analysis_functions.R")
out <- "./DESeq2"
if(!dir.exists(out)) {dir.create(out)}


#################### Import and Format ########################


### For testing

#path <- "/Volumes/GoogleDrive/My Drive/Projects/Tanya_Beluga_RNAseq_R008586/"
meta <- read_excel(file.path(path,"RNAseq_beluga_comparison_Group_description.xlsx"),1)
meta <- meta[,1:3]
meta <- mutate(meta,Sub_group_SLE=str_replace(Sub_group_SLE,"-","Not_SLE"))
### add column for batch
batch1 <- c(
"DLB-1521","DLB-1536",
"DLB-1615","DLB-1626",
"DLB-1653")

meta <- mutate(meta,batch=ifelse(Name %in% batch1,"batch1","batch2"))


dplyr::count(meta,Group)
# Group  n
# 1 Arctic 34
# 2    SLE 39
# 3    Tag 12

dplyr::count(meta,Sub_group_SLE)

# Sub_group_SLE  n
# 1        Center  5
# 2    Downstream 16
# 3       Not_SLE 46
# 4      Upstream 18

meta <- column_to_rownames(meta,"Name")

## from genpipes
counts <- read.delim(file.path(path,"rawCountMatrix.csv"),check.names = F)
annotation <- counts[,c(1,2)]

rownames(counts) <- counts$Gene
counts <- counts[,-c(1,2)]

## get gene descriptions
## dataset in format of hsapiens_gene_ensembl
ensembl <- useEnsembl(biomart = "genes",dataset = "dleucas_gene_ensembl",version = 105)
descriptions <- getBM(attributes = c("ensembl_gene_id","description"),filters = "ensembl_gene_id",rownames(counts),mart = ensembl)
descriptions <- descriptions[match(rownames(counts),descriptions$ensembl_gene_id),]
annotation$description <- descriptions$description

counts <- counts[,match(rownames(meta),colnames(counts))]
# colnames(counts)[!colnames(counts) %in% rownames(meta)]
# rownames(meta)[!rownames(meta) %in% colnames(counts)]
# ## change TagHI-19-06_JA" to TagHI-19-06_JA
# colnames(counts)[colnames(counts) == "TagHI-19-06_JA"] <- "TagHI-19-06"
# counts <- counts[,match(rownames(meta),colnames(counts))]
# 

comparisons <- data.frame(group2=c("SLE","Tag"),group1=rep("Arctic",2))


 contrasts <- map(1:nrow(comparisons), function(x) {
   treat <- comparisons[x,1]
   control <- comparisons[x,2]

   return(c("Group",treat,control))
})

names(contrasts) <- paste(comparisons[,1],comparisons[,2],sep = "_vs_")


################ Differential Expression analysis #####################

################## Full data set ######################################
dds <- DESeqDataSetFromMatrix(colData = meta,
                              countData = counts,
                              design= ~ Group)
mcols(dds) <- cbind(mcols(dds), annotation)

##filtering
total <- dim(dds)[1]
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
kept <- dim(dds)[1]
filtered = total - kept
# 349 filtered
## 21646 kept



### Visualization ############
dds <- DESeq(dds)

vsd <- vst(dds)
p <- plotPCA2(vsd,intgroup = "Group")
ggsave(p,filename = file.path(out,"PCA.png"))




out <- "."
sample_heatmap(vsd,"Group",file.path(out,"Sample_heatmap2.png"))

top_genes_heatmap(vsd, top = 50,annotation_cols = "Group",file=file.path(out,"top_50_most_variable_genes_heatmap.png"))


DESeq2::plotDispEsts(dds)

normd_counts <- counts(dds,normalized = T)
normd_counts <- rownames_to_column(as.data.frame(normd_counts),"GeneID")
#write.csv(normd_counts, file = file.path(out,"Muscle_Stem_cells_normalized_counts.csv"),row.names=F)


############# Testing #####################
##this will do all comparisons
Res <- lapply(contrasts,function(x) {
  res <- results(dds,contrast = x,alpha = 0.05)
  lfcShrink(dds=dds,res = res,type = "ashr")
})


map2(Res,names(Res),~maPlot.lists(.x,.y,alpha = 0.05))




Res_with_annot <- lapply(Res,function(x) {
  
  dat <- as.data.frame(x)
  dat <- rownames_to_column(dat,"Ensembl_Gene_ID")
  descriptions <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),filters = "ensembl_gene_id",values=dat$Ensembl_Gene_ID,mart = ensembl)
  dat <- merge(dat,descriptions,by.x="Ensembl_Gene_ID",by.y="ensembl_gene_id")
  colnames(dat)[8] <- "Gene_Symbol" 
  return(dat)
})


file_names <- file.path(out,paste0(names(Res_with_annot),"_DEtable.csv"))
map2(Res_with_annot,file_names,function(x,y) {
  write.csv(x,file = y,row.names=F)
})

Sig_genes <- map(Res_with_annot,~filter(.x,!is.na(padj),padj<=0.05,abs(log2FoldChange)>=1))




######### Gene Ontology #############################

library(gprofiler2)
gene_lists <- map(Sig_genes,"Ensembl_Gene_ID")

map2(Sig_genes,names(Sig_genes), ~write.table(.x,file=file.path(out,paste0(.y,"FC2_fdr_0.05_sig_gene_names.txt")),
                                              row.names=F,col.names=F,quote=F))

### use portal
# go_results <- map(gene_lists, ~gost(query = .x, 
#                                     organism = "dleucas", ordered_query = F, 
#                                     multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
#                                     measure_underrepresentation = FALSE, evcodes = TRUE, 
#                                     user_threshold = 0.05, correction_method = "g_SCS", 
#                                     domain_scope = "known", custom_bg = NULL, 
#                                     numeric_ns = "", sources = c("GO","KEGG")))
# 
# 


## format results
gp <- read.csv("gProfiler_dleucas_SLE_vs_Artic_2-23-2022_4-29-45 PM__intersections.csv")


###################### Sub Group SLE Analysis #####################


## Sub group
```{r}
dds_SLE <- dds[,colData(dds)$Group == "SLE"]
colData(dds_SLE)$Sub_group_SLE <- as.factor(colData(dds_SLE)$Sub_group_SLE)
design(dds_SLE) <- ~Sub_group_SLE
dds_SLE <- estimateSizeFactors(dds_SLE)

### Visualization #########################################
dds_SLE <- DESeq(dds_SLE)

vsd_SLE <- vst(dds_SLE)

p <- plotPCA2(vsd_SLE,intgroup = "Sub_group_SLE")
ggsave(p,filename = file.path(out,"PCA_Sub_group_SLE.png"),width=7,height = 2)


sample_heatmap(vsd_SLE,"Sub_group_SLE",file.path(out,"Sample_heatmap_Sub_group_SLE.png"))



top_genes_heatmap(vsd_SLE, top = 50, annotation_cols="Sub_group_SLE",file=file.path(out,"top_50_most_variable_genes_heatmap_Sub_group_SLE.png"))


DESeq2::plotDispEsts(dds_SLE)



comparisons <- data.frame(group2=c("Upstream","Downstream","Upstream"),group1=c(rep("Center",2),"Downstream"))

### Testing ################################
contrasts <- map(1:nrow(comparisons), function(x) {
  treat <- comparisons[x,1]
  control <- comparisons[x,2]
  
  return(c("Sub_group_SLE",treat,control))
})


names(contrasts) <- paste(comparisons[,1],comparisons[,2],sep = "_vs_")

##this will do all comparisons
Res <- lapply(contrasts,function(x) {
  res <- results(dds_SLE,contrast = x,alpha = 0.05)
  lfcShrink(dds=dds_SLE, res=res,type = "ashr")
})


map2(Res,names(Res),~maPlot.lists(.x,.y,alpha = 0.05))



Res_with_annot <- lapply(Res,function(x) {
  
  dat <- as.data.frame(x)
  dat <- rownames_to_column(dat,"Ensembl_Gene_ID")
  descriptions <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),filters = "ensembl_gene_id",values=dat$Ensembl_Gene_ID,mart = ensembl)
  dat <- merge(dat,descriptions,by.x="Ensembl_Gene_ID",by.y="ensembl_gene_id")
  colnames(dat)[8] <- "Gene_Symbol" 
  return(dat)
})

file_names <- file.path(out,paste0(names(Res_with_annot),"_DEtable.csv"))
map2(Res_with_annot,file_names,function(x,y) {
  write.csv(x,file = y,row.names=F)
})

Sig_genes <- map(Res_with_annot,~filter(.x,!is.na(padj),padj<=0.05,abs(log2FoldChange)>=1))

map2(Sig_genes,names(Sig_genes), ~write.table(.x,file=file.path(out,paste0(.y,"FC2_fdr_0.05_sig_gene_names.txt")),
                                              row.names=F,col.names=F,quote=F))


#### Visualize #############
custom_genes_heatmap(vsd_SLE,annotation_cols="Group",gene_list = fcntrl_vs_Treatment.topgenes,file= "test.pdf",
                     col_font=10,samples=NULL,
                     show_rownames=F,show_colnames,row_font=5,treeheight_row = 0,scale="row",cellwidth=10)


#################### Save and export #####################################

saveRDS(dds,file=file.path(out,"full_DESeq2_DataSet.rds"))
saveRDS(dds_SLE,file=file.path(out,"SLE_DESeq2_DataSet.rds"))

saveRDS(vsd,file=file.path(out,"vsd.rds"))
saveRDS(vsd_SLE,file=file.path(out,"SLE_vsd.rds"))

capture.output(sessionInfo(),file = file.path(out,"SessionInfo.txt"))



