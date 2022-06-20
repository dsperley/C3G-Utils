library(edgeR)
library(ggplot2)
library(reshape2)
library(gplots)
library(data.table)

mytheme <-  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
theme_set(mytheme)

#setwd("/Volumes/GoogleDrive/My Drive/Projects/Tanya_Beluga_RNAseq_R008586/")
source("./RNA_DE_analysis_functions.R")
out <- "./Limma_Voom"
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
meta <- mutate(meta,across(everything(),as.factor))

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


design <- model.matrix(~ 0 + Group, data = meta)
colnames(design) <- gsub("Group","",colnames(design))
#is.fullrank(design)


### cant take a data frame of contrasts
comparisons <- data.frame(group2=c("SLE","Tag"),group1=rep("Arctic",2))
contrasts <- map_chr(1:nrow(comparisons), function(i) {paste0(comparisons[i,1],"-",comparisons[i,2])})
contr.matrix <- makeContrasts(contrasts=x,levels =colnames(design))


################ Differential Expression analysis #####################
dge <-  DGEList(counts, samples = meta, group = meta$Group , genes = annotation)
total <- nrow(dge$counts)
#[1] 21995

keep <- rowSums(dge$counts)>=10
sum(keep)
#[1] 21646
dge <- dge[keep,]
dge <- calcNormFactors(dge)
v <- voom(dge, design , plot = FALSE)

### QC EDA #########
##### library sizes #######
col.cell <-
  c("#d73027",
    "#fc8d59",
    "#fee090",
    "#33a02c",
    "#91bfdb",
    "#4575b4")[meta$Group]
png(file=file.path(out,"Library_sizes.png"),width = 10,height = 8,units = 'in',res = 300)
barplot(
  dge$samples$lib.size,
  #names = anno$Condition,
  las = 2,
  col = col.cell
)
dev.off()

##### PCA ###############
p <- plotPCA2(v,intgroup="Group")
ggsave(p,filename = file.path(out,"PCA.png"),width = 10,height = 8,units = 'in',dpi = 300)

##### MDS ###############
p <- plotMDS2(dge,intgroup = "Group")
ggsave(p,filename = file.path(out,"MDS.png"),width = 10,height = 8,units = 'in',dpi = 300)


##### BCV ###############
dge <- estimateDisp(dge, design)
png(file = file.path(out,"Dispersion_estimate.png"),width = 10,height = 8,units = 'in',res = 300)
plotBCV(dge, cex = 0.5, main = "Dispersion Estimate")
dev.off()


##### sample heatmaps ######
sample_heatmap(v,"Group",file=file.path(out,"sample_heatmap.pdf"))

### Testing #################
fit <- lmFit(v, design)
plotSA(fit, main = "voom: Mean-variance trend")

fit2 <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit2)

png(file = file.path(out, "Mean-variance_trend.png"),width = 10,height = 8,units = 'in',res = 300)
plotSA(efit, main = "Final model: Mean-variance trend")
dev.off()

#######Extracting results ##########

################Summary ###################
summary.table <- as.data.frame.matrix(summary(decideTests(efit, 
                                                          p.value = aPVALcutoff, lfc = log2FCcutval)))
summary.table <- t(summary.table)
newdf <-
  rbind(summary.table, data.frame(Down = paste0("*FC >= ", FCcutoff, "*FDR >= ", aPVALcutoff),NotSig="",Up=""))

newdf <- rownames_to_column(newdf,"Group")
newdf[newdf$Group == "1","Group"] <- ""
write.csv(newdf, file.path(out, "table_summary_DEG_stat.csv"),
          row.names = F)
png(
  file = file.path(out, "table_summary_DEG_stat.png"),
  width = 10, height = 8,units = 'in',res = 300)

grid.table(newdf)
dev.off()

############# DE tables ########################

contrast_list <- colnames(contr.matrix)

my_contrast_list <- list()




all_res <- map(contrast_list,function(c) {
  ##tcnt <- tcnt + 1
  res <- data.frame(topTable(efit,coef = contrast_list[y],
    adjust = "BH",n = nrow(dge$counts)))
   res <- dplyr::select(res,-AveExpr,-t,-B)
   write.csv(res,file.path(out,paste0(c,"_DEG_TotalValidGenes_in_Genome.csv")),row.names=F)
   return(res)
})
names(all_res) <- colnames(contr.matrix)

deg_res <- map(contrast_list,function(c) {
    ##tcnt <- tcnt + 1
    res <- data.frame(topTable(efit,coef = contrast_list[y], p.value = aPVALcutoff,
                               adjust = "BH",n = nrow(dge$counts,lfc = log2FCcutval,sort.by="p")))
    res <- dplyr::select(res,-AveExpr,-t,-B)
    write.csv(res,file.path(out,paste0(c,"_DiffExpressed_FDR_",aPVALcutoff,"_FC_",FCcutoff,".csv")),row.names=F)
    return(res)

})
  
################# Visualization of Results ################

### DEG heatmaps #####
map2(deg_res,file,function(x,f) {
gene_list <- if(nrow(deg)<50) {
  rownames(deg)
} else {
  rownames(deg)[1:50]
}  
custom_genes_heatmap(dge,annotation_cols="Group",gene_list = gene_list,file= f,
                     col_font=10,samples=NULL,
                     show_rownames=F,show_colnames,row_font=5,treeheight_row = 0,scale="row",cellwidth=10)
})

##### Violin Plots #################

###### Venn
  
