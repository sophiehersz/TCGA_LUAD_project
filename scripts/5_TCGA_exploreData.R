# Exploration and visualization of RNAseq data

setwd("/Users/sophie/Dropbox (The Francis Crick)/Sophie/Astra_project/Repositories/TCGA_LUAD_project/")

library(readr)
library(dplyr)
library(gridExtra)
library(umap)
library(ggplot2)
library(Rtsne)

# Functions ####################################################################
dimReduction <- function(normalized_counts, sample_info, method="UMAP"){
  #' Returns dataframe with projection coordinates and accompanying sample info
  #' @param normalized_counts: dataframe containing normalized RNAseq counts for all samples
  #' @param sample_info:  dataframe with the sample barcodes and their corresponding mutation_status and sample_type
  #' @param method: string specifying the projection method to be used (can take "UMAP" or "TSNE" values)
  
  print("Scaling counts...")
  # transpose counts_table to have genes as columns and patients as rows
  normalized_counts <- t(normalized_counts)
  # apply Z-score standardization to counts
  scaled_counts <- scale(normalized_counts)
  
  # dimensionality reduction
  if (method == "UMAP"){
    print("Computing UMAP...")
    set.seed(42)
    projection <- umap(scaled_counts)
    print("Done")
    # convert results into dataframe (in preparation for plotting)
    projection_df <- data.frame(projection$layout)
    
  } else if (method == "TSNE"){
    print("Computing tSNE...")
    set.seed(42) 
    projection <- Rtsne(scaled_counts)
    print("Done")
    projection_df <- data.frame(projection$Y)
  }
  
  # Add info to daraframe
  projection_df <- projection_df %>% 
    mutate(mutation_status = as.factor(sample_info$mutation_status),
           tumor = as.factor(sample_info$tumor),
           cases = row.names(.),
           plate =  substr(cases, 22, 25))
  
  return(projection_df)
}   

plotProjection <- function(data, title=NULL){
  #' Plots the projection results
  #' @param data: dataframe returned by dimReduction() function

  ggplot(data, aes(x=X1, y=X2, color = mutation_status)) + geom_jitter() +
    ggtitle(title) +
    scale_color_manual(values=c("#2b83ba", "#d7191c"), "STK11 mutation") +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 15),
          axis.title = element_blank(), 
          axis.text = element_text(size = 12),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"))
}

# Load data ####################################################################
DE_results <- read_csv("results/DE_analysis/ExactTest_pc_genes_tumor_only_no_outliers.csv")
TMM_counts <- read_csv('data/data_for_ML/luad_rna_TMM_counts_pc_genes_tumor_only_no_outliers.csv')
sample_info <- read_csv('data/data_for_ML/sample_info_tumor_only_no_outliers.csv')


# Plot results tables ##########################################################
# filter genes by FDR and logFC
de_genes <- DE_results %>% 
  filter(FDR <= 0.01) %>% 
  filter(logFC > 1 | logFC < -1) %>% 
  arrange("FDR")
# round values for visualization purposes
de_genes[,3:4] <- round(de_genes[,3:4], digits = 3)
de_genes[,5:6] <- signif(de_genes[,5:6], digits=3)
# make table with top20 genes for visualization, arranged by FDR
dev.off()
grid.table(de_genes[1:20,], rows=NULL)


# check significance of genes of interest
kaufman_genes <- c('AVPI1', 'BAG1', 'CPS1', 'DUSP4', 'FGA', 'GLCE', 'HAL', 'IRS2', 
                   'MUC5AC', 'PDE4D', 'PTP4A1', 'RFK', 'SIK1', 'TACC2', 'TESC', 'TFF1')
k_genes <- DE_results %>% 
  filter(Symbol %in% kaufman_genes)
k_genes[,3:4] <- round(k_genes[,3:4], digits = 3)
k_genes[,5:6] <- signif(k_genes[,5:6], digits=3)
dev.off()
grid.table(k_genes, rows=NULL)


# Dimensionality reduction ####################################################
umap_df <- dimReduction(TMM_counts, sample_info, method="UMAP")
plotProjection(umap_df, title="UMAP projection of LUAD RNAseq samples")
    
tsne_df <- dimReduction(TMM_counts, sample_info, method="TSNE")
plotProjection(tsne_df, title="t-SNE projection of LUAD RNAseq samples")

