# Prepares RNAseq counts_table for DE analysis from TCGA RNAseq HTseq counts data
# If required, also filters the table to keep only genes of interest (e.g protein coding)

setwd("/Users/sophie/Dropbox (The Francis Crick)/Sophie/Astra_project/Repositories/TCGA_LUAD_project")

library(readr)
library(biomaRt)

# Functions ####################################################################
getSampleInfo <- function(patient_list, rna_df){
  #' Returns dataframe with barcodes and corresponding mutation status and sample type (tumor vs normal) for all RNAseq samples
  #' belonging to the patient IDs in patient_list
  #' Also reorders dataframe to place WT samples first
  #' @param patient_list a dataframe containing patient IDs "cases.submitter_id" and corresponding "mutation_status"
  #' @param rna_df a dataframe from GDC containing patient IDs "cases.submitter_id", "sample_type" and sample barcodes "cases"
  
  # add mutation status to all samples and keep only samples belonging to patient_list
  df <- full_join(rna_df[, c("cases.submitter_id", "sample_type", "cases", "id")], 
                  patient_list, by = "cases.submitter_id") %>% 
    filter(cases.submitter_id %in% patient_list$cases.submitter_id)
  
  # recodes sample_type as 0 or 1 (normal or tumor, respectively)
  # also sets the mutation_status of all normal tissue samples as 0
  # finally, drops redundant columns and reorders rows to place WT samples first
  df <- df %>% 
    mutate(tumor = case_when(sample_type == 'Solid Tissue Normal' ~ 0, TRUE ~ 1)) %>% 
    mutate(mutation_status = case_when(sample_type == 'Solid Tissue Normal' ~ 0, TRUE ~ mutation_status)) %>% 
    relocate(cases) %>% 
    dplyr::select(-sample_type) %>% 
    arrange(mutation_status)
  
  return(df)
}

getCountsTable <- function(path, sample_info, rna_df){
  
  #' Prepares and returns a counts_table dataframe for DE analysis
  #' @param path path to the GDCdata folder where counts files are stored
  #' @param sample_info a dataframe with the names of the counts files to be considered ("cases", e.g. "TCGA-44-7660-01A-11R-2066-07")
  #' @param rna_df a dataframe from GDC containing patient IDs ("case.submitter_id") and names of counts files for all patients ("cases")
  
  # initialize counts_table
  counts_table <- data.frame(ensembl_gene_id = character())
  
  # list of folder names and sample barcodes
  folders <-sample_info$id
  sample_barcodes <- sample_info$cases
  
  # fill counts_table with data
  for (i in seq_along(folders)) {
    file <- list.files(paste(path, "/", folders[i], sep = ""), pattern = "*.gz") # get file name
    sample <- sample_barcodes[i] # get sample barcode
    
    data <- read.delim(paste(path, "/", folders[i], "/", file, sep = ""), 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       header=FALSE, 
                       col.names=c("ensembl_gene_id", sample))
    
    counts_table <- full_join(counts_table, data, by = "ensembl_gene_id") # adds data to counts_table
    
    print(paste("Fetching RNAseq data... ", i, "/", length(folders)))
  }
  
  # counts_table cleanup
  rows_to_delete <- which(substring(counts_table$ensembl_gene_id, 1, 1) != "E") # Remove non-gene rows (rows not named after an ensembl gene ID ENSG)
  counts_table <- counts_table[-rows_to_delete,]
  gene_ids <- sub("\\..*", "", counts_table$ensembl_gene_id) # remove versions from ensembl gene IDs
  counts_table$ensembl_gene_id <- gene_ids
  
  print('Job done')
  
  return(counts_table)
}

filterGenes <- function(counts_table, mart, gene_types=NULL){
  #' Adds gene symbols to counts table and selects genes of interest
  #' @param counts_table: dataframe with RNAseq counts for all patients, created with getCountsTable function
  #' Needs to contain "ensembl_gene_id" column and other columns are named after sample barcodes
  #' @param mart: Mart object containing gene symbols and type (retrieved wih biomaRt package)
  #' @param gene_types: an array of strings with the gene_biotypes to be considered (e.g. "protein_coding"). If gene_types
  #' is specified, counts_table is filtered to keep only genes belonging to gene_types. If not, only gene symbols will be
  #' added to counts_table (no filtering will be performed)
  
  # Use Mart object to get info for genes in counts_table, and add it to the table
  genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype', 'description'),
                 filters = 'ensembl_gene_id', 
                 values = counts_table$ensembl_gene_id,
                 mart = mart)
  
  # adds gene symbols and biotypes to counts_table
  counts_table_genes <- left_join(counts_table, genes[, c('ensembl_gene_id', 'external_gene_name', "gene_biotype")], 
                                  by = 'ensembl_gene_id') 
  
  # filters genes if requested
  if(!is.null(gene_types)){
    counts_table_genes <- counts_table_genes %>% 
      filter(gene_biotype %in% gene_types) %>% 
      relocate(external_gene_name) %>% 
      dplyr::select(-gene_biotype)
  }
  
  else {
    counts_table_genes <- counts_table_genes %>% 
      relocate(external_gene_name) %>%
      dplyr::select(-gene_biotype)
  }
  
  return(counts_table_genes)
  
}


# Load data ####################################################################
path_to_rna_data = "data/GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/" # path where RNAseq counts data is stored
patient_list <- read_csv("data/luad_patient_list.csv") # list of patient IDs to be considered
luad_rna <- read_csv("data/luad_rna.csv") # dataframe containing the file names of the RNAseq data and corresponding patient IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") # database with gene information

# creates sample_info dataframe, with info on mutation_status and sample_type, for all RNAseq samples belonging to patient_list
sample_info <- getSampleInfo(patient_list, rna_df=luad_rna)
write_csv(sample_info, "data/luad_rna_sample_info.csv")

# Prepare counts_table #########################################################
# creates a "general" dataframe containing RNAseq counts for all genes and all samples in sample_info
counts_table <- getCountsTable(path=path_to_rna_data, sample_info, rna_df=luad_rna)
write_csv(counts_table, "data/luad_raw_counts_table.csv")

# filters counts_table to keep only the genes of interest
counts_table_pc <- filterGenes(counts_table, mart, gene_types=c("protein_coding"))
write_csv(counts_table_pc, "data/luad_raw_counts_table_pc.csv")
