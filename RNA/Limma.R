# Objective:
#  to reanalyze to 2 microarray studies and compare the gene expression profiles between lession and non-lesional skin tissues, lesional and heatlhy control skin tissue
# 
# Datasets:
#   GSE57178: 
#     wheal/lesional and non-lesional skin samples from Patients 1, 2, 3, 4 6, 7 and the healthy control samples Normal 11-15.
#   GSE72542: wheal/lesional and non-lesional skin samples from patients A7, A8, A15, A17, A18, A19 (the non-lesional samples for this patient is incorrectly labelled in GEO as healthy control), A20, A22, XX and the healthy control samples from B1, B2, B6, B7, B8, B9, B10. Additionally, 
#           the blood samples from patients A1-5, A7, A8, A10, A12-A22 and AXX and control blood from B1-10.


#library(BiocManager)
library(magrittr)
library(dplyr)
library(GEOquery)
library(ggplot2)
#library(data.table)
library(oligo)
library(limma)
library(org.Hs.eg.db)
library(stringr)

#setwd("~/Documents/C3G/Projects/Chronic_Spontaneous_Urticaria_GEO_Project_R004021")

## functions
toFold  <- function(x){
  x[x >= 0] = 2^(x[x >= 0])
  x[x < 0] = -2^(-x[x < 0])
  return(x)
}


plotdf <- function(DE_list,groups,ref) {
  ## Takes a list of DE results from both studies, extracts the expression values, 
  ## creates a data frame suitable for plotting
  map(DE_l_vs_nl,~dplyr::select(.x,ENSEMBL,starts_with("GSM"))) %>% 
    purrr::reduce(function(x,y) {inner_join(x,y,by="ENSEMBL")}) %>%
    pivot_longer(-ENSEMBL,names_to = "Sample",values_to = "Expression") %>%
    extract(Sample,into=c("Study","Group"),regex="(GSM[0-9]+).(.*)",remove=F) %>%
    mutate(Group=str_replace_all(Group,c("normal_skin"="healthy","lesional_skin" = "lesional")),
           Study=str_replace_all(Study,c("GSM13.*"= "GSE57178","GSM18.*"="GSE72542"))) %>%
    filter(Group %in% groups) %>%
    mutate(Group=factor(Group) %>% relevel(ref = ref)) %>%
    group_by(ENSEMBL,Study) %>%
    mutate(Z_expression=scale(Expression))
}


genePlot <- function(df,g,meta) {
  
  gene_symbol <- goi[goi$ensembl_gene_id == g,"external_gene_name"]
  ens <- g
  FC <- meta[ meta$ENSEMBL == g,"metafc"] %>% toFold()
  p <- meta[ meta$ENSEMBL == g,"metap"] 
  padj <- meta[ meta$ENSEMBL == g,"FDR"]
  
  
  title <- sprintf("%s (%s)\n FC: %.2f  P.value: %.2e  adj.P.value: %.2e",gene_symbol,ens,FC,p,padj)
  
  p <- filter(df,ENSEMBL == g) %>%
    ggplot(aes(Group,Z_expression,color=Group)) +
    geom_boxplot(aes(color=Group)) +
    geom_jitter(aes(color=Group),width = 0.1) +
    labs(y="Expression (Z-score)",title=title) +
    theme_bw()
  
  return(p)
}




##############
##  GSE57178
##############
gse = "GSE57178"
eset = getGEO(gse,destdir = "./")[[1]]
pData(eset)$"SampleID" = sampleNames(eset)
pData(eset)$"Group" = pData(eset)$"Tissue" =  pData(eset)$"tissue:ch1" %>% gsub(" ","_",.) %>% gsub("-","_",.)
pData(eset)$"Patient" = pData(eset)$"patient:ch1"
eset <- eset[,pData(eset)$"Patient" != "8"]
fData(eset)$"ENSEMBL" = fData(eset)$"mrna_assignment" %>% str_extract_all("gene:ENSG[0-9]*") %>%    sapply(unique,simplify=F) %>% sapply(paste,collapse=" /// ")  %>% gsub("gene:","",.)

plotMDS(eset,pch=19,col=as.numeric(as.factor(pData(eset)$Tissue)))
legend("topleft",pch=19,legend = levels(as.factor(pData(eset)$Tissue)), 
                          col=palette()[1:nlevels(as.factor(pData(eset)$Tissue))])

pdf(file= "GSE57178_MDS.pdf")
plotMDS(eset,pch=19,col=as.numeric(as.factor(pData(eset)$Tissue)))
legend("topleft",pch=19,legend = levels(as.factor(pData(eset)$Tissue)), 
       col=palette()[1:nlevels(as.factor(pData(eset)$Tissue))])
dev.off()

#DE
Group = pData(eset)$"Group"
design = model.matrix(~-1+Group); colnames(design) %<>% gsub("^Group","",.)
coefs = c("lesional_skin-non_lesional_skin","lesional_skin-normal_skin")
cont.matrix = makeContrasts(contrasts=coefs, levels=design)
fit = lmFit(eset,design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

lapply(coefs,function(c) {
  top = topTable(fit2,coef=c,number=Inf,genelist=NULL,confint = T)
  top$"Rank" = 1:nrow(top)
  top$"FC" = toFold(top$"logFC")
  top <- cbind(fData(eset), top[featureNames(eset),])
  x = exprs(eset); colnames(x) %<>% paste(pData(eset)$Group,sep=":")
  top %<>% cbind(x)
  write.csv(top,file = paste0(c,"_GSE57178_DE.csv"))
})
  


################


#A7, A8, A15, A17, A18, A19 (the non-lesional samples for this patient is incorrectly labelled in GEO as healthy control), A20, A22, XX and the healthy control samples from B1, B2, B6, B7, B8, B9, B10. Additionally, 
#           the blood samples from patients A1-5, A7, A8, A10, A12-A22 and AXX and control blood from B1-10.

gse = "GSE72542" #
eset = getGEO(gse,destdir = "./")[[1]]
pData(eset)$"SampleID" = sampleNames(eset)
pData(eset)$"Group" = pData(eset)$"Tissue" =  pData(eset)$"tissue:ch1" %>% 
                      gsub("wheal","lesional",.) %>%
                      gsub("^skin$","healthy",.) %>%
                      gsub("-","_",.) %>%
                      gsub(" ","_",.)
pData(eset)$"Patient" = pData(eset)$"individual id:ch1"
##
pData(eset)[pData(eset)$Patient == "A19"& pData(eset)$Tissue == "healthy",c("Tissue","Group")] <- rep("non_lesional",2)

Patients_keep <- c("A7", "A8", "A15", "A17", "A18", "A19","A20", "A22", "XX","B1", "B2", "B6", "B7", "B8", "B9", "B10")
tissues_rm <- "blood"
eset <- eset[,pData(eset)$"Patient" %in% Patients_keep & pData(eset)$Tissue != tissues_rm]
## another fix
pData(eset)[pData(eset)$Tissue == "non_lesional",c("Group","Tissue")] <- rep("non_lesional_skin",2)

library(biomaRt)
ensembl<- useEnsembl(biomart = 'genes', 
                     dataset = 'hsapiens_gene_ensembl',
                     version = "GRCh37")

filter <- c("agilent_sureprint_g3_ge_8x60k_v2")
attributes <-c("agilent_sureprint_g3_ge_8x60k_v2","ensembl_gene_id")
values <- fData(eset)$SPOT_ID

bm <- getBM(attributes = attributes,filters = filter,mart = ensembl,values=values)
## deal with probes assigning to multiple Ensembl gene ids
bm_dedup <-bm$ensembl_gene_id %>% 
  split(bm$agilent_sureprint_g3_ge_8x60k_v2) %>%
  map(unique) %>% sapply(paste,collapse=";") %>% as.data.frame()

fData(eset)$ENSEMBL <- bm_dedup[match(fData(eset)$SPOT_ID,rownames(bm_dedup)),] 




plotMDS(eset,pch=19,col=as.numeric(as.factor(pData(eset)$Tissue)))
legend("bottomleft",pch=19,legend = levels(as.factor(pData(eset)$Tissue)), 
       col=palette()[1:nlevels(as.factor(pData(eset)$Tissue))])



pdf(file= "GSE72542_MDS.pdf")
plotMDS(eset,pch=19,col=as.numeric(as.factor(pData(eset)$Tissue)))
legend("topleft",pch=19,legend = levels(as.factor(pData(eset)$Tissue)), 
       col=palette()[1:nlevels(as.factor(pData(eset)$Tissue))])
dev.off()

#DE
Group = pData(eset)$"Group"
design = model.matrix(~-1+Group); colnames(design) %<>% gsub("^Group","",.)
coefs = c("lesional-non_lesional_skin","lesional-healthy")
cont.matrix = makeContrasts(contrasts=coefs, levels=design)
fit = lmFit(eset,design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

lapply(coefs,function(c) {
  top = topTable(fit2,coef=c,number=Inf,genelist=NULL,confint = T)
  top$"Rank" = 1:nrow(top)
  top$"FC" = toFold(top$"logFC")
  top <- cbind(fData(eset), top[featureNames(eset),])
  x = exprs(eset); colnames(x) %<>% paste(pData(eset)$Group,sep=":")
  top %<>% cbind(x)
  write.csv(top,file = paste0(c,"_GSE72542_DE.csv"))
})

#########################################################
################# Meta Analysis #########################
##########################################################

### Curation, making the gene lists match
library(MetaVolcanoR)
DE_lesional_vs_non_lesional_files <- list.files(pattern = "lesional.*non_lesional_skin_GSE.*DE.csv")
DE_l_vs_nl <- map(DE_lesional_vs_non_lesional_files,read.csv)
names(DE_l_vs_nl) <- c("Affy","Agilent")


## select probes/probsets assigned only to one gene, then select probe/probeset most differentially expressed per gene
DE_l_vs_nl<- map(DE_l_vs_nl, ~filter(.x,str_count(ENSEMBL,"ENS")==1))
DE_l_vs_nl <- map(DE_l_vs_nl, ~group_by(.x,ENSEMBL) %>%
  dplyr::slice(which.min(adj.P.Val)))


DE_l_vs_nl$Agilent <- dplyr::rename(DE_l_vs_nl$Agilent,TX_ID = ENSEMBL_ID)

common_genes <- intersect(DE_l_vs_nl$Affy$ENSEMBL,DE_l_vs_nl$Agilent$ENSEMBL)
DE_l_vs_nl <- map(DE_l_vs_nl,~filter(.x,ENSEMBL %in% common_genes))


### Meta
meta_l_vs_nl_comb <- combining_mv(DE_l_vs_nl,pcriteria = "P.Value",
                             foldchangecol = "LogFC",
                             metafc = "Mean",geneidcol = "ENSEMBL")

meta_l_vs_nl_comb_res <- meta_l_vs_nl_comb@metaresult
meta_l_vs_nl_comb_res$FDR <- p.adjust(meta_l_vs_nl_comb_res$metap,method = "BH")

filter <- c("ensembl_gene_id")
attributes <-c("ensembl_gene_id","external_gene_name","description")
values <- meta_l_vs_nl_comb_res$ENSEMBL

bm <- getBM(attributes = attributes,filters = filter,mart = ensembl,values=values)
meta_l_vs_nl_comb_res <- merge(meta_l_vs_nl_comb_res,bm,by.x="ENSEMBL",by.y="ensembl_gene_id")
meta_l_vs_nl_comb_res$idx <- NULL

meta_l_vs_nl_comb_res <- map(DE_l_vs_nl,~dplyr::select(.x,ENSEMBL,starts_with("GSM"))) %>% 
  purrr::reduce(function(x,y) {inner_join(x,y,by="ENSEMBL")}) %>%
  inner_join(meta_l_vs_nl_comb_res,.,by="ENSEMBL") 

write.csv(meta_l_vs_nl_comb_res,file = "lesional_vs_non_lessional_meta_DE_results.csv",row.names = F)



#########################################

DE_lesional_vs_healthy_files <- c("lesional_skin-normal_skin_GSE57178_DE.csv","lesional-healthy_GSE72542_DE.csv")
DE_l_vs_h <- map(DE_lesional_vs_healthy_files,read.csv)
names(DE_l_vs_h) <- c("Affy","Agilent")


## select probes/probsets assigned only to one gene, then select probe/probeset most differentially expressed per gene
DE_l_vs_h<- map(DE_l_vs_h, ~filter(.x,str_count(ENSEMBL,"ENS")==1))
DE_l_vs_h <- map(DE_l_vs_h, ~group_by(.x,ENSEMBL) %>%
                    dplyr::slice(which.min(P.Value)))


DE_l_vs_h$Agilent <- dplyr::rename(DE_l_vs_h$Agilent,TX_ID = ENSEMBL_ID)

common_genes <- intersect(DE_l_vs_h$Affy$ENSEMBL,DE_l_vs_h$Agilent$ENSEMBL)
DE_l_vs_h <- map(DE_l_vs_h,~filter(.x,ENSEMBL %in% common_genes))


### Meta
meta_l_vs_h_comb <- combining_mv(DE_l_vs_h,pcriteria = "P.Value",
                                  foldchangecol = "LogFC",
                                  metafc = "Mean",geneidcol = "ENSEMBL")


meta_l_vs_h_comb_res <- meta_l_vs_h_comb@metaresult
meta_l_vs_h_comb_res$FDR <- p.adjust(meta_l_vs_h_comb_res$metap,method = "BH")

filter <- c("ensembl_gene_id")
attributes <-c("ensembl_gene_id","external_gene_name","description")
values <- meta_l_vs_h_comb_res$ENSEMBL

bm <- getBM(attributes = attributes,filters = filter,mart = ensembl,values=values)
meta_l_vs_h_comb_res <- merge(meta_l_vs_h_comb_res,bm,by.x="ENSEMBL",by.y="ensembl_gene_id")
meta_l_vs_h_comb_res$idx <- NULL

### add expression values

meta_l_vs_h_comb_res <- map(DE_l_vs_h,~dplyr::select(.x,ENSEMBL,starts_with("GSM"))) %>% 
  purrr::reduce(function(x,y) {inner_join(x,y,by="ENSEMBL")}) %>%
  inner_join(meta_l_vs_h_comb_res,.,by="ENSEMBL") 
write.csv(meta_l_vs_h_comb_res,file = "lesional_vs_healthy_meta_DE_results.csv",row.names = F)



#############################################################
## genes of interest
genes <- c("IL4", "IL4R", "IL13", "IL5", "GATA3", "STAT6", "CCR8", "IL25",
"IL17A", "IL17F", "IL21", "IL22", "IL23A", "SOCS3", "STAT3", "IL6", "IL12A", "RORC", "IL21R",
"IFNG", "TBX21", "IL2", "TNF", "STAT4",
"IL10", "FOXP3","STAT6", "TGFB1", "IL1R1", "STAT6", "EBI3",
"PDGFRA", "IL33", "IL9", "IL18")


filter <- c("external_gene_name")
attributes <-c("ensembl_gene_id","external_gene_name","description")
values <- genes

bm <- getBM(attributes = attributes,filters = filter,mart = ensembl,values=values)
## are these in dataset?
goi <- bm[bm$ensembl_gene_id %in% DE_l_vs_h$Affy$ENSEMBL,]
not <- genes[!genes %in% goi$external_gene_name] 
## "RORC"   "TNF"    "FOXP3"  "PDGFRA"


##################################################### 
########### Plot



l_vs_nl_df <- plotdf(DE_l_vs_nl,c("lesional","non_lesional"),"non_lesional")


pdf(file="lesional_vs_non_lesional_Z_scores.pdf")
map(goi$ensembl_gene_id,~genePlot(l_vs_nl_df,.x,meta_l_vs_nl_comb_res))
dev.off()

l_vs_h_df <- plotdf(DE_l_vs_h,c("lesional","healthy"),"healthy")


pdf(file="lesional_vs_healthy_Z_scores.pdf")
map(goi$ensembl_gene_id,~genePlot(l_vs_h_df,.x,meta_l_vs_h_comb_res))
dev.off()


l_vs_h_Z <- ungroup(l_vs_h_df) %>%
dplyr::select(ENSEMBL,Sample,Z_expression) %>% 
pivot_wider(names_from = Sample,values_from = Z_expression) 

write.csv(l_vs_h_Z,file = "lesional_vs_healthy_Z_scores.csv",row.names = F)

l_vs_nl_Z <- ungroup(l_vs_nl_df) %>%
  dplyr::select(ENSEMBL,Sample,Z_expression) %>% 
  pivot_wider(names_from = Sample,values_from = Z_expression) 

write.csv(l_vs_nl_Z,file = "lesional_vs_non_lesional_Z_scores.csv",row.names=F)


