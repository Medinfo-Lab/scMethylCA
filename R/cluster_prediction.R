#' RNA standard data
#'
#' @param file_tmp_length The length of each gene entered
#' @param method Methods for data standardization, TPM, FPKM, CPM
#' @param file_tmp Enter the RNA expression matrix
#'
#' @return Expression matrix after normalization
#' @export
#'
#' @examples
RNA_data_standard <- function(file_tmp,file_tmp_length,method="tpm"){
  file_data <- STAD_seurat_data
  file_gene_length <- STAD_expr_length_choose

  if (method=="tpm") {
    file_data_count <- file_data
    file_gene_length_kb <- file_gene_length[,1]/1000

    rpk <- file_data_count/file_gene_length_kb
    tpm_matrix <- t(t(rpk) / colSums(rpk) * 1e6)
    return(tpm_matrix)
  }else if(method=="fpkm"){
    file_data_count <- file_data
    file_gene_length_kb <- file_gene_length[,1]/1000

    rpk <- file_data_count/file_gene_length_kb
    fpkm_matrix <- t(t(rpk) / colSums(rpk) * 1e9)
    return(fpkm_matrix)
  }else if(method=="cpm"){
    file_data_count <- file_data
    total_counts <- colSums(file_data_count)

    cpm_matrix <- t(t(file_data_count) / total_counts) * 1e6
    return(cpm_matrix)
  }
}


#' LASSO gene
#'
#' @param file_tmp Enter the RNA expression matrix
#' @param file_group_tmp Grouping information for the sample
#' @param test_size Training-to-testing ratio
#'
#' @return Returns the expression matrix after selection
#' @export
#'
#' @examples
LASSO_feature <- function(file_tmp,file_group_tmp,test_size=0.3){
  source_python_file <- system.file(package = "scMethylCA")
  source_def <- paste(source_python_file,"python/lasso_feature.py",sep = "/")
  reticulate::source_python(source_def)
  # reticulate::source_python("inst/python/lasso_feature.py")

  result <- reticulate::py$lasso_feature(file_tmp,file_group_tmp,test_size)

  return(result)
}



#' Grouping model prediction
#'
#' @param file_tmp Enter the RNA expression matrix
#' @param test_size Training-to-testing ratio
#' @param model Model selection
#' @param file_group Grouping information for the sample
#'
#' @return Post-training data
#' @export
#'
#' @examples
Group_model_predict <- function(file_tmp,file_group,test_size=0.3,model="lgbm"){
  source_python_file <- system.file(package = "scMethylCA")
  source_def <- paste(source_python_file,"python/model_algorithm.py",sep = "/")
  reticulate::source_python(source_def)
  # reticulate::source_python("inst/python/model_algorithm.py")

  result <- reticulate::py$model_algorithm(file_tmp,file_group,test_size,model)

  return(result)
}
