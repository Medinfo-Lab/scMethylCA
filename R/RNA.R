#' UMAP dimensionality reduction data was used for quasi-temporal analysis
#'
#' @param rds_tmp RDS objection
#'
#' @return monocle3 object and Pseudotime data
#' @export
#'
#' @examples
Pseudotime_analysis <- function(rds_tmp){
  pbmc <- rds_tmp
  expression_matrix <- GetAssayData(pbmc,assay="RNA",slot="counts")
  cell_metadata <- pbmc@meta.data
  gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix),
                              row.names = rownames(expression_matrix))
  cds <- new_cell_data_set(expression_data = expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)
  #preprocess_cds==seurat NormalizeData+ScaleData+RunPCA
  # cds <- preprocess_cds(cds)
  # cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP")

  int.embed <- Embeddings(pbmc,reduction = "umap")
  # int.embed <- int.embed[rownames(cds.embed),]
  reducedDims(cds)$UMAP <- int.embed


  cds <- cluster_cells(cds,reduction_method = "UMAP")
  cds <- learn_graph(cds,
                     learn_graph_control = list(
                       minimal_branch_len = 1,
                       euclidean_distance_ratio = 1))
  cds <- order_cells(cds)

  pseudotime_values <- pseudotime(cds, reduction_method = "UMAP")
  pseudo_seurat <- as.data.frame(pseudotime_values)
  pseudo_seurat_UMAP <- cds@int_colData$reducedDims$UMAP
  pseudo_seurat_UMAP <- as.data.frame(pseudo_seurat_UMAP)
  pseudo_seurat$UMAP1 <- pseudo_seurat_UMAP$umap_1
  pseudo_seurat$UMAP2 <- pseudo_seurat_UMAP$umap_2
  pseudo_seurat$group <- pbmc$seurat_clusters

  pseudo_list <- list()
  pseudo_list$pseudo_object <- cds
  pseudo_list$pseudo_data <- pseudo_seurat
  return(pseudo_list)
}


#' Differential analysis of RNA expression matrices
#'
#' @param file_tmp RNA expression matrix
#' @param file_group Group Information
#' @param group_id1 Group Character1
#' @param group_id2 Group Character2
#' @param method_model Statistical methods
#'
#' @return Difference analysis results
#' @export
#'
#' @examples
RNA_group_variance_analysis <- function(file_data,file_group,group_id1,group_id2=NULL,method_model="wilcoxon"){
  # file_data <- file_tmp
  # file_group <- TCGA_STAD_normalized_group

  if (is.null(group_id2)) {
    group_id_sample <- subset(file_group, group == group_id1)
    group_noid_sample <- subset(file_group, group != group_id1)
  }else{
    group_id_sample <- subset(file_group, group == group_id1)
    group_noid_sample <- subset(file_group, group == group_id2)
  }

  meth_level_group_id_sample <- file_data[,group_id_sample$sample]
  meth_level_group_noid_sample <- file_data[,group_noid_sample$sample]

  meth_level_group_id_sample <- as.matrix(meth_level_group_id_sample)
  meth_level_group_noid_sample <- as.matrix(meth_level_group_noid_sample)
  file_data_rowMean <- rowMeans(file_data,na.rm = T)


  if (method_model == "wilcoxon") {
    wilcoxon_data <- data.frame(row.names = rownames(file_data))
    wilcoxon_data$AveExpr <- file_data_rowMean

    for (i in 1:nrow(meth_level_group_id_sample)) {
      tmp <- wilcox.test(meth_level_group_id_sample[i,],meth_level_group_noid_sample[i,])
      wilcoxon_data$logFC[i] <- log2(mean(meth_level_group_id_sample[i,],na.rm = T)/mean(meth_level_group_noid_sample[i,],na.rm = T))
      wilcoxon_data$P.value[i] <- tmp$p.value
      wilcoxon_data$adj_P.value[i] <- p.adjust(tmp$p.value,method = "BH")
    }
    return(wilcoxon_data)
  }else if(method_model == "var"){
    var_data <- data.frame(row.names = rownames(file_data))
    var_data$AveExpr <- file_data_rowMean

    for (i in 1:nrow(meth_level_group_id_sample)) {
      tmp <- var.test(meth_level_group_id_sample[i,],meth_level_group_noid_sample[i,])
      var_data$logFC[i] <- log2(mean(meth_level_group_id_sample[i,],na.rm = T)/mean(meth_level_group_noid_sample[i,],na.rm = T))
      var_data$p_value[i] <- tmp$p.value
      var_data$adj_P.value[i] <- p.adjust(tmp$p.value,method = "BH")
    }
    return(var_data)
  }else if(method_model == "t.test"){
    t.test_data <- data.frame(row.names = rownames(file_data))
    t.test_data$AveExpr <- file_data_rowMean

    for (i in 1:nrow(meth_level_group_id_sample)) {
      tmp <- t.test(meth_level_group_id_sample[i,],meth_level_group_noid_sample[i,])
      t.test_data$logFC[i] <- log2(mean(meth_level_group_id_sample[i,],na.rm = T)/mean(meth_level_group_noid_sample[i,],na.rm = T))
      t.test_data$p_value[i] <- tmp$p.value
      t.test_data$adj_P.value[i] <- p.adjust(tmp$p.value,method = "BH")
    }
    return(t.test_data)
  }
}


#' Enrichment analysis of RNA expression matrices
#'
#' @param file_tmp Expression matrix derived from differential analysis
#' @param method Enrichment methods
#'
#' @return Enrichment results
#' @export
#'
#' @examples
RNA_enrichment_analysis <- function(file_tmp,method="GO"){
  file_data <- FF_DEG_Dx_choose

  entrezid_all = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
                        keys = rownames(file_data), #将输入的gene_name列进行数据转换
                        keytype = "SYMBOL", #输入数据的类型
                        column = "ENTREZID")#输出数据的类型

  GO_database <- 'org.Hs.eg.db'
  KEGG_database <- 'hsa'

  if(method=="GO"){
    GO<-enrichGO(entrezid_all$ENTREZID,#GO富集分析
                 OrgDb = GO_database,
                 keyType = "ENTREZID",#设定读取的gene ID类型
                 ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                 pvalueCutoff = 1,#设定p值阈值
                 qvalueCutoff = 1,#设定q值阈值
                 readable = T)
    return(GO)
  }else if(method=="KEGG"){
    KEGG<-enrichKEGG(entrezid_all$ENTREZID,#KEGG富集分析
                     keyType = "kegg",
                     organism = KEGG_database,
                     pAdjustMethod = 'fdr',  #指定p值校正方法
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
    return(KEGG)
  }
}
