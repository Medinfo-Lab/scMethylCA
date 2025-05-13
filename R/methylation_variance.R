#' distinguish chrXX:xxx-xxx as chrXX, xxx, xxx
#'
#' @param file_tmp Chromosome physical segment data, chrXX:xxx-xxx format or chrXX, xxx, xxx format
#'
#' @return Processed physical fragment data, chrXX, xxx, xxx format or chrXX:xxx-xxx format
#' @export
#'
#' @examples
Chr_region_process <- function(file_tmp,method){
  if(method=="split"){
    chr_data <- file_tmp
    chr_data_frame <- data.frame(region=chr_data)
    chr_split <- chr_data_frame %>%
      separate(region, into = c("chr", "start", "end"), sep = ":|-")
    chr_split$start <- as.numeric(chr_split$start)
    chr_split$end <- as.numeric(chr_split$end)
    return(chr_split)
  }else if(method=="paste"){
    chr_data <- file_tmp
    chr_data_paste <- sprintf("%s:%d-%d",chr_data$chr,chr_data$start,chr_data$end)
    chr_data_paste_frame <- as.data.frame(chr_data_paste)
    colnames(chr_data_paste_frame) <- "chr"
    return(chr_data_paste_frame)
  }
}





#' Chromosome physical segment annotation
#'
#' @param file_tmp Chromosome physical segment data, chrXX, xxx, xxx format
#' @param txdb_file Genome annotation databases, such as TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' @return Annotated chromosome physical segment data
#' @export
#'
#' @examples
Chr_region_annot <- function(file_tmp,txdb_file){
  peaks_data <- file_tmp

  peaks <- GRanges(
    seqnames = peaks_data$chr,
    ranges = IRanges(start = peaks_data$start, end = peaks_data$end)
    # strand = Rle(strand = c("+", "-", "+"))
  )
  txdb <- txdb_file
  peakAnno <- annotatePeak(peaks, TxDb = txdb, annoDb="org.Hs.eg.db", addFlankGeneInfo=T)
  peakAnno_frame <- as.data.frame(peakAnno)
  return(peakAnno_frame)
}






#' Group analysis
#'
#' @param file_tmp DNA methylation or chromatin accessibility level data
#' @param file_group Group Information
#' @param group_id1 Group Character1
#' @param group_id2 Group Character2
#' @param mothed Statistical methods
#'
#' @return Difference analysis results
#' @export
#'
#' @examples
Meth_group_variance_analysis <- function(file_data,file_group,group_id1,group_id2=NULL,method="t.test"){
  # file_data <- STAD_CpG_methlevel_choose
  # file_group <- STAD_sample_choose_methlevel
  # group_id1 <- "Gastric Cancer Line"

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


  if (method=="wilcoxon") {
    wilcoxon_data <- data.frame(row.names = rownames(file_data))
    wilcoxon_data$AveMethLevel <- file_data_rowMean

    for (i in 1:nrow(meth_level_group_id_sample)) {
      tmp <- wilcox.test(meth_level_group_id_sample[i,],meth_level_group_noid_sample[i,])
      wilcoxon_data$logFC[i] <- log2(mean(meth_level_group_id_sample[i,],na.rm = T)/mean(meth_level_group_noid_sample[i,],na.rm = T))
      wilcoxon_data$P_value[i] <- tmp$p.value
      wilcoxon_data$adj_P.value[i] <- p.adjust(tmp$p.value,method = "BH")
    }
    return(wilcoxon_data)
  }else if(method=="t.test"){
    t_test_data <- data.frame(row.names = rownames(file_data))
    t_test_data$AveMethLevel <- file_data_rowMean

    for (i in 1:nrow(meth_level_group_id_sample)) {
      tmp <- t.test(meth_level_group_id_sample[i,],meth_level_group_noid_sample[i,])
      t_test_data$logFC[i] <- log2(mean(meth_level_group_id_sample[i,],na.rm = T)/mean(meth_level_group_noid_sample[i,],na.rm = T))
      t_test_data$P_value[i] <- tmp$p.value
      t_test_data$adj_P.value[i] <- p.adjust(tmp$p.value,method = "BH")
    }
    return(t_test_data)
  }else if(method=="var"){
    var_data <- data.frame(row.names = rownames(file_data))
    var_data$AveMethLevel <- file_data_rowMean

    for (i in 1:nrow(meth_level_group_id_sample)) {
      tmp <- var.test(meth_level_group_id_sample[i,],meth_level_group_noid_sample[i,])
      var_data$logFC[i] <- log2(mean(meth_level_group_id_sample[i,],na.rm = T)/mean(meth_level_group_noid_sample[i,],na.rm = T))
      var_data$P_value[i] <- tmp$p.value
      var_data$adj_P.value[i] <- p.adjust(tmp$p.value,method = "BH")
    }
    return(var_data)
  }
}



#' GSEA analysis
#'
#' @param file_tmp Difference analysis results
#' @param gmt_data user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
#' @param pvaluecutoff adjusted pvalue cutoff
#'
#' @return gseaResult object
#' @export
#'
#' @examples
GSEA_analysis <- function(file_tmp,gmt_file,pvaluecutoff=1){
  file_data <- file_tmp
  gmt_data=read.gmt(gmt_file)

  file_data_trim <- file_data[,c("seqnames","start","end","ENSEMBL","SYMBOL","logFC","P_value","adj_P.value")]
  rt=file_data_trim[order(file_data_trim[,"logFC"],decreasing=T),]
  logFC=as.vector(rt[,"logFC"])
  names(logFC)=as.vector(rt[,"SYMBOL"])
  GSEA_data=GSEA(logFC,TERM2GENE = gmt_data,pvalueCutoff = pvaluecutoff)
  return(GSEA_data)
}
