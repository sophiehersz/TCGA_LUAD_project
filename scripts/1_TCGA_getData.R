# Queries the GDC database to retrieve RNAseq (raw counts) and mutation data for a given project
# Mutation data from the MC3 Pan-Cancer study (https://gdc.cancer.gov/about-data/publications/mc3-2017)

library(readr)
library(dplyr)
library(TCGAbiolinks)
library(TCGAmutations)

setwd("/Users/sophie/Dropbox (The Francis Crick)/Sophie/Astra_project/Repositories/TCGA_LUAD_project/")

# Functions ####################################################################  
getGeneExp <- function(project, download = FALSE, destination_folder) {
  #' Uses the TCGAbiolinks package (https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) to query the GDC database.
  #' Returns a dataframe with a list of RNAseq raw counts files available (non-normalized, HTseq counts) for the given project
  #' Also has the option to download the data - in that case, destination_folder must be provided
  #' @param project: a string corresponding to the name of the TCGA project of interest (e.g. 'TCGA-LUAD'). Multiple entries accepted.
  #' @param download: boolean, whether or not to download the data from the query
  #' @param destination_folder: if download=TRUE, the destination folder where the data will be saved
  
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    legacy = FALSE
  )
  
  exp_df <- getResults(query)
  
  if (download) {
    setwd(destination_folder)
    GDCdownload(query_exp, method = "api")
  }
  
  return(exp_df)
}

getMut <- function(project) {
  #' Gets MAF file from the MC3 Pan-Cancer study (https://gdc.cancer.gov/about-data/publications/mc3-2017)
  #' Uses the package TCGAmutations (https://github.com/PoisonAlien/TCGAmutations)
  #' Returns a dataframe with somatic mutation data for the given project
  #' @param project: a string corresponding to the name of the project of interest (e.g. 'TCGA-LUAD')
  
  maf <- tcga_load(study = project)
  maf_df <- maf@data %>% 
    mutate(cases.submitter_id = substr(Tumor_Sample_Barcode, 1, 12)) # adds a column with patient IDs
  
  return(maf_df)
}

# Get and save RNAseq data #####################################################
luad_rna <- getGeneExp(project = "TCGA-LUAD")
write_csv(luad_rna, file = "data/luad_rna.csv")

# Get and save mutation data ###################################################
luad_mut <- getMut("LUAD")
write_csv(luad_mut, file = "data/luad_mut_mc3.csv")


