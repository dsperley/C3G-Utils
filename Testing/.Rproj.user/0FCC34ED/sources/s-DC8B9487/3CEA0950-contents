


############################################
## add project specific functions here


##############################################

## This Function will allow us to make a MA plot from a list of results objects
maPlot.lists <- function(x,i,alpha=0.05) {
  pdf(file=file.path(out,paste0(i, '_maPlot.pdf', sep='')))
  plotMA(x, main=paste(i, 'alpha = ',alpha), 
         alpha=alpha, ylim=c(-6,6))
  abline(h=c(1,-1), col='red')
  dev.off()
}

# ##This function will add either count means for fpkm means to a DE table
# addMeans<-function(x,means,comp) {
#   detable<-as.data.frame(x)
#   
#   ##get the two group names
#   groups<-unlist(strsplit(comp,"vs"))
#   group_means<-means[,groups]
#   
#   ##match the gene names in the mean df with the results df
#   group_means<-group_means[match(rownames(detable),rownames(group_means)),]
#   detable<-cbind(detable,group_means)
#   return(detable)}


## modified DESeq2 plotPCA function, to plot additional pc dimensions
plotPCA2 <- function (object, intgroup = "condition", ntop = 500, dim1= 1,dim2 = 2,returnData = FALSE,addLabels=F,...) 
{
  ## ... is geom_text params if addLabels is True
  #### deal with different objects
  if (class(object) == "DESeqTransform") {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t(assay(object)[select,])
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  
  sampleNames <- colnames(object)
  
  }
  
  if (class(object) == "EList") {
    rv <- rowVars(object$E)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    mat <- t(object$E[select,])  
  
    if (!all(intgroup %in% names(object$targets))) {
      stop("the argument 'intgroup' should specify columns of y$samples")
    }
    
    intgroup.df <- object$targets[,intgroup, drop = FALSE]
    sampleNames <- colnames(object$E) 
  }
  
 ######
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    intgroup.df[[intgroup]]
  }
  
  
  pca <- prcomp(mat)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  d <- data.frame( pca$x[, dim1],pca$x[, dim2], group = group, 
                   intgroup.df, name = sampleNames)
  dim_names <- paste0("PC",c(dim1,dim2))
  colnames(d)[1:2] <- dim_names
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(dim1,dim2)]
    return(d)
  }
  p <- ggplot(data = d, aes_string(x = dim_names[1], y = dim_names[2], color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0(dim_names[1],": ", round(percentVar[dim1] * 100), "% variance")) + 
    ylab(paste0(dim_names[2],": ", round(percentVar[dim2] * 100), "% variance"))
 
  if(addLabels == T) {
    p <- p +
      geom_text(...)
  }
  
  return(p)
    
}


plotMDS2 <- function(object,intgroup,addLabels=F,...) {
  if (class(object) == "EList") {
    
    if (!all(intgroup %in% names(object$targets))) {
      stop("the argument 'intgroup' should specify columns of y$samples")
    }
    mds <- plotMDS(object)
    plot_df  <-  data.frame(
      Dim1 = mds$x,
      Dim2 = mds$y,
      RepID = rownames(object$targets),
      Group = object$targets[,intgroup])
    
    } else if (class(object) == "DGEList") {
      if (!all(intgroup %in% names(object$samples))) {
        stop("the argument 'intgroup' should specify columns of y$samples")
      }  
      mds <- plotMDS(object)
      plot_df  <-  data.frame(
                Dim1 = mds$x,
                Dim2 = mds$y,
                RepID = rownames(object$samples),
                Group = object$samples[,intgroup])
  
    }
  
  p <- ggplot(plot_df, aes(Dim1, Dim2, colour = Group)) +
    geom_point(stroke = 1.5, size = 3,alpha = 1) +
    xlab("Leading logFC dim 1") +
    ylab("Leading logFC dim 2") +
    if(addLabels == T) {
      p <- p +
        geom_text(...)
    }
   return(p)
        
  }    


sample_heatmap <- function(object,annotation_cols=NULL,file) {
  
  if (class(object) == "DESeqTransform") {
    if(!is.null(annotation_cols)) {
      if(!annotation_cols %in% colnames(colData(object))) {
        stop("annotation variable is not in colData")
      }
    }
    mat <- t(assay(object))
    anno_df <- colData(object)[,annotation_cols,drop=F]
    anno_df <- as.data.frame(anno_df)
    
  } else if (class(object) == "EList") {
    if(!is.null(annotation_cols)) {
      if(!annotation_cols %in% colnames(object$targets)) {
        stop("annotation variable is not in targets data frame")
      }
    }
    mat <- t(object$E)
    anno_df <- object$targets[,annotation_cols,drop=F]
  
  }
    
  dist_mat <- as.matrix(dist(mat))
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  if (is.null(annotation_cols)) {
    pheatmap(dist_mat,color = colors,cellwidth=10, cellheight=10)
  } else {
    pheatmap(dist_mat,color = colors,annotation_col = anno_df,filename = file,cellwidth=10, cellheight=10)
  }
}

custom_genes_heatmap <- function(object,meta,gene_list=NULL,annotation_cols=NULL,samples=NULL,file,...) {
  require("pheatmap")
  ## this function will take a list of gene ids - can be genes of interest, DEG's ect and
  ## generate a heatmap of those gene expression values
  ## for DESeq2 obje is rld or vst
  if("DESeqDataObject" %in% class(object)) {
  
    if(!is.null(annotation_cols)) {
      if(!annotation_cols %in% colnames(colData(object))) {
        stop("annotation variable is not in colData")
      }
    }
  
    if(is.null(samples)) {
     ## use all samples
      samples <- colnames(object)
    } else if (sum(samples %in% colnames(object) == 0)) {
      stop("samples not found in object")
    }
      
    
    
    if(is.null(gene_list)) {
      stop("Must provide a gene list")
    } else if(sum(gene_list %in% rownames(assay(object)))!=length(gene_list)) {
      warning("Not all genes in list")
    }else {  
    select <-  which(gene_list %in% assay(object))  
    anno_df <- colData(object)[samples,annotation_cols,drop=F]
    anno_df <- as.data.frame(anno_df)
    mat <- assay(vsd)[select,samples]
    }
    
  } else if ("ExpressionSet" %in% class(object)) {
    
      if(!is.null(annotation_cols)) {
        if(!annotation_cols %in% colnames(pData(object))) {
          stop("annotation variable is not in colData")
        }
      }
    
    if(is.null(samples)) {
      ## use all samples
      samples <- colnames(object)
    } else if (sum(samples %in% colnames(object) == 0)) {
      stop("samples not found in object")
    }
    
  
    if(is.null(gene_list)) {
      stop("Must provide a gene list")
    } else if(sum(gene_list %in% featureNames(object))!=length(gene_list)) {
      warning("Not all genes in list")
    }else {  
       
      anno_df <- pData(object)[samples,annotation_cols,drop=F]
      anno_df <- as.data.frame(anno_df)
      mat <- exprs(object)[gene_list,samples]
      
    }
  }  
  else if("matrix" %in% class(object)) {
    if(!is.null(annotation_cols)) {
      if(!annotation_cols %in% colnames(meta)) {
        stop("annotation variable is not in col")
      }
    }
    
    if(is.null(samples)) {
      ## use all samples
      samples <- colnames(object)
    } else if (sum(samples %in% colnames(object) == 0)) {
      stop("samples not found in object")
    }
    
    
    if(is.null(gene_list)) {
      stop("Must provide a gene list")
    } else if(sum(gene_list %in% rownames(object))!=length(gene_list)) {
      warning("Not all genes in list")
    }else {  
      
      anno_df <- meta[match(samples,rownames(meta)),annotation_cols,drop=F]
      anno_df <- as.data.frame(anno_df)
      mat <- object[gene_list,match(samples,colnames(object))]
      
    }
  }
    
  else if("DGEList" %in% class(object)) {
    
    if(!is.null(annotation_cols)) {
      if(!annotation_cols %in% colnames(object$samples)) {
        stop("annotation variable is not in samples")
      }
    }
    
    if(is.null(samples)) {
      ## use all samples
      samples <- colnames(object)
    } else if (sum(samples %in% colnames(object) == 0)) {
      stop("samples not found in object")
    }
    
    
    if(is.null(gene_list)) {
      stop("Must provide a gene list")
    } else if(sum(gene_list %in% rownames(object))!=length(gene_list)) {
      warning("Not all genes in list")
    }else {  
      
      anno_df <- object$samples[match(samples,rownames(object$samples)),annotation_cols,drop=F]
      lcpms <- cpm(object,log=T)
      mat <- lcpms[gene_list,match(samples,colnames(lcpms))]
      
    }
  }
  
  
  
  if (is.null(annotation_cols)) {
    pheatmap(mat, filename = file,...)
  } else {
    pheatmap(mat, annotation_col = anno_df,filename=file,...) 
  }
}


Violin_plot <- function(gene,meta=NULL,object,group=NULL,DE=NULL) {

  
  if("ExpressionSet" %in%  class(object)) {
  values <- scale(log2(exprs(object)[gene,] + 1))  
  plot_df <- data.frame(sample=rownames(values),expression=values)
  ## add annotation
  plot_df <- merge(plot_df,pData(object),by.x="sample",by.y="title")
    if(!is.null(DE)) {
      res <- DE[DE$Gene == gene,]
      gene_symbol <- gene
      ens <- res["ensembl_gene_id"] 
      FC <- res["FC"]
      p <- res["P.Value"]
      padj <- res["adj.P.Val"]
      
      title <- sprintf("%s (%s)\n FC: %.2f  P.value: %.2e  adj.P.value: %.2e",gene_symbol,ens,FC,p,padj)
    }
  } else if("matrix" %in%  class(object)) {
    values <- scale(log2(object)[gene,] + 1) 
    plot_df <- data.frame(sample=rownames(values),expression=values)
    ## add annotation
    plot_df <- merge(plot_df,meta,by.x="sample",by.y="Title")
    if(!is.null(DE)) {
      res <- DE[DE$Gene == gene,]
      gene_symbol <- gene
      ens <- res["ensembl_gene_id"] 
      FC <- res["FC"]
      p <- res["P.Value"]
      padj <- res["adj.P.Val"]
      
      title <- sprintf("%s (%s)\n FC: %.2f  P.value: %.2e  adj.P.value: %.2e",gene_symbol,ens,FC,p,padj)
    }   
  }
  
  
    
   
  g <- ggplot(plot_df,aes_string(group,"expression")) +
    geom_violin(aes_string(fill=group)) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_boxplot(aes_string(fill=group),width=0.05) +
    ylab(paste(gene,"expression Z-score")) +
    guides(fill=FALSE) +
    #stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05 )
    theme_classic()

   if(!is.null(DE)) {
     g <- g + ggtitle(title)
   }
  
  return(g)

  
}

top_genes_heatmap <- function(object,top=50,annotation_cols=NULL,row_font=5,file) {
  
  if (class(object) == "DESeqTransform") {
    if(!is.null(annotation_cols)) {
      if(!annotation_cols %in% colnames(colData(object))) {
        stop("annotation variable is not in colData")
      }
    }
    select <- order(rowVars(assay(object)),decreasing=TRUE)[1:top]
    mat <- t(assay(object)[select,])
    anno_df <- colData(object)[,annotation_cols,drop=F]
    anno_df <- as.data.frame(anno_df)
    
  }
  
  
  if (is.null(annotation_cols)) {
    pheatmap(mat,cellwidth=10, cellheight=10, show_rownames=T,scale = "row",treeheight_row = 0,show_colnames = F,border_color = NA,fontsize_row = row_font,filename = file)
  } else {
    pheatmap(mat, cellwidth=10, cellheight=10,show_rownames=T,scale = "row",treeheight_row = 0,show_colnames = F,border_color = NA,annotation_col = anno_df,fontsize_row =row_font,filename = file)
  }
}

##### pull out top most variable genes
get_top_genes <- function(object,top=50) {
  if (class(object) == "DESeqTransform") {
    select <- order(rowVars(assay(object)),decreasing=TRUE)[1:top]
    assay(object)[select,]
    }
  
  
}


## DE is Res_with_anno, gene lengths is a dataframe with 2 columns,- GeneID, and gene length, genome -"mm10", and "ensGene"
run_goseq <- function(de,gene_lengths,genome,ids) {
  assayed.genes <- de$Ensembl_Gene_ID
  de.genes <- assayed.genes[de$padj<=0.05]
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  
  bias.data <- gene_lengths$V2
  names(bias.data) <- gene_lengths$V1
  bias.data <- bias.data[match(names(gene.vector),names(bias.data))]
  pwf <- nullp(gene.vector,genome,ids,bias.data = bias.data)
  go <- goseq(pwf,genome,ids)
  return(go)
}

