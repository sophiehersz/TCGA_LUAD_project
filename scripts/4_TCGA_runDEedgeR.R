# Performs edgeR DE analysis

setwd("/Users/sophie/Dropbox (The Francis Crick)/Sophie/Astra_project/Repositories/TCGA_LUAD_project/")

library(readr)
library(dplyr)
library(edgeR)
library(tidyr)
library(stringr)

# Functions ####################################################################
filterData <- function(counts_table, sample_info, tumor_only=TRUE, outliers=NULL){
  #' filters counts_table to exclude normal samples and/or outlier samples if necessary
  #' @param counts_table: a dataframe with RNAseq raw counts, created with getCountsTable function. 
  #' @param sample_info: a dataframe with the sample barcodes and their corresponding sample_type
  #' @param tumor_only: booleanm whether or not to exclude normal samples
  #' @param outliers: dataframe with outlier sample barcodes to be excluded, if any
  
  # get list of normal sample barcodes
  if (tumor_only){
    normal_samples <- sample_info %>% 
      filter(tumor == 0) %>% 
      pull(cases)
  } else {
    normal_samples = NULL}
  
  # get list of outlier sample barcodes
  if(!is.null(outliers)){
    outlier_samples <- outliers %>% 
      pull(cases)
  } else {
    outlier_samples = NULL
  }
  
  # columns to be deleted from counts_table
  to_delete <- which(colnames(counts_table) %in% c(normal_samples, outlier_samples))
  
  # delete samples to be excluded from counts_table and also from sample_info
  counts_table_filtered <- counts_table[,-to_delete]
  sample_info_filtered <- sample_info %>% 
    filter(cases %in% colnames(counts_table_filtered))
  
  # check if sample order is identical between counts_table and sample_info
  counts_samples <- colnames(counts_table_filtered[,3:ncol(counts_table_filtered)])
  info_samples <- sample_info_filtered$cases
  test = identical(counts_samples, info_samples)
  
  print(paste("Identical tables: ", toString(test)))
  
  assign("counts_table_filtered", counts_table_filtered, .GlobalEnv)
  assign("sample_info_filtered", sample_info_filtered, .GlobalEnv)
  
}

prepareDGEList <- function(counts_table, sample_info, norm_method="TMM", filterByExp = TRUE){
  #' prepares DGEList object for DE analysis
  #' also performs count normalization for subsequent dimensionality reduction
  #' @param counts_table: a dataframe with RNAseq raw counts, created with getCountsTable function. 
  #' Sample barcodes as columns. Columns 1 and 2 correspond to "gene_symbols" and "ensembl_gene_ids" respectively
  #' @param sample_info: a dataframe with the sample barcodes and their corresponding mutation_status and sample_type
  #' @param norm_method: a string specifying the normalization method to be used (for CalcNormFactors {edgeR} function)
  #' @param filterByExp: boolean, specifies whether or not to run the default filterByExp {edgeR} function, to discard genes
  #' with low counts (keep genes with CPM ≥ 1 in at least 70% of samples per group)
  
  groups <- sample_info$mutation_status
  
  # create DGEList object
  y <- DGEList(counts = counts_table[,3:ncol(counts_table)], genes = counts_table[,2], group = groups)
  y$genes$Symbol <- counts_table$external_gene_name
  
  # keep gene with CPM ≥ 1 in at least 70% of samples (per group)
  if (filterByExp) {
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    dim(y)
  }

  # Remove duplicated gene symbols by keeping the entries with the highest counts
  o <- order(rowSums(y$counts), decreasing=TRUE)
  y <- y[o,]
  d <- duplicated(y$genes$Symbol)
  y <- y[!d,]
  
  # Re-compute library sizes after the filtering
  y$samples$lib.size <- colSums(y$counts)
  
  # Rename y$counts and y$genes rows as Ensembl IDs and remove the genes$genes column
  rownames(y$counts) <- rownames(y$genes) <- y$genes$genes
  y$genes$genes <- NULL
  
  # Normalize counts
  y <- calcNormFactors(y, method = norm_method)
  
  return(y)
  
}

DEAnalysis <- function(y){
  #' performs edgeR DE analysis using the exact test method (comparison between 2 groups)
  #' returns a dataframe containing logFC, p-values and FDR for all genes
  #' @param y a DGEList object containing the group variable
  
  # creates design matrix
  design <- model.matrix(~0+group, data=y$samples)
  colnames(design) <- c("WT","MUT")
  
  # estimates dispersion (to account for technical and biological variability)
  print("Estimating dispersion...")
  y <- estimateDisp(y, design)
  
  # performs exact test
  print("Performing exact test...")
  et <- exactTest(y, pair=c(1,2))
  assign("et",et,.GlobalEnv)
  
  # generates table containing all genes and their respective test values, ordered by FDR
  all_genes <- topTags(et, n = dim(y)[1])$table
  
  return(all_genes)
}

# load data ####################################################################
counts_table_pc <- read_csv("data/luad_raw_counts_table_pc.csv") # protein-coding genes only
sample_info <- read_csv("data/luad_rna_sample_info.csv")
outliers <- read_csv("data/luad_rna_outliers.csv") # outlier samples to be removed from dataset, if any

# Filter data ##################################################################
# option to remove normal tissue samples and/or outlier samples
filterData(counts_table_pc, sample_info, tumor_only=TRUE, outliers=outliers)
write_csv(sample_info_filtered, 'data/data_for_ML/sample_info_tumor_only_no_outliers.csv')
write_csv(counts_table_filtered, 'data/data_for_ML/luad_raw_counts_table_tumor_only_no_outliers.csv')

# run DE analysis ##############################################################
y <- prepareDGEList(counts_table_filtered, sample_info_filtered)
all_genes_et <- DEAnalysis(y)
write_csv(all_genes_et, "results/DE_analysis/ExactTest_pc_genes_tumor_only_no_outliers.csv")

# get normalized counts table for further analysis
TMM_counts <- cpm(y, normalized.lib.sizes=TRUE) # uses normalized library sizes, so that TMM normalized counts values are returned
TMM_counts <- data.frame(TMM_counts, check.names = FALSE)
write_csv(TMM_counts, 'data/data_for_ML/luad_rna_TMM_counts_pc_genes_tumor_only_no_outliers.csv')

# get gene names (order preserved with rows in TMM_counts dataframe)
genes <- data.frame(Symbol = y$genes$Symbol)
write_csv(genes, 'data/data_for_ML/genes.csv')


# OPTIONAL #####################################################################
# run analysis including  normal samples #######################################
filterData(counts_table_pc, sample_info, tumor_only=FALSE, outliers=outliers)
write_csv(sample_info_filtered, 'data/data_for_ML/sample_info_tumor_normal_no_outliers.csv')

y <- prepareDGEList(counts_table_filtered, sample_info_filtered)
all_genes_et <- DEAnalysis(y)
write_csv(all_genes_et, "results/DE_analysis/ExactTest_pc_genes_tumor_normal_no_outliers.csv")

TMM_counts <- cpm(y, normalized.lib.sizes=TRUE) # uses normalized library sizes, so that TMM normalized counts values are returned
TMM_counts <- data.frame(TMM_counts, check.names = FALSE)
write_csv(TMM_counts, 'data/data_for_ML/luad_rna_TMM_counts_pc_genes_tumor_normal_no_outliers.csv')

# # save data for ML
# genes <- data.frame(Symbol = y$genes$Symbol)
# write_csv(genes, 'data/sophie_ML/genes.csv')



