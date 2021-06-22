# TCGA_LUAD_project
Work in progress
## Data sources

### RNAseq data - HTseq raw counts
Obtained from the GDC Data Portal (https://portal.gdc.cancer.gov/).

Downloaded using the `TCGAbiolinks` R package (https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html).

### Mutation data
Obtained from the MC3 Pan-cancer study (https://gdc.cancer.gov/about-data/publications/mc3-2017).

Downloaded using the `TCGAmutations` R package (https://github.com/PoisonAlien/TCGAmutations)

### WSI data
Image embeddings downloaded from https://github.com/binli123/dsmil-wsi

Publication: *Li et al., 2021* (https://arxiv.org/abs/2011.08939)

---
## RNAseq analysis pipeline
Differential expression (DE) analysis performed using edgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html)

![alt text](https://github.com/sophiehersz/TCGA_LUAD_project/blob/main/thumbnails/RNAseq_pipeline.png?raw=true)

## Running the pipeline
### Data download
For RNAseq and mutation data download, run `1_TCGA_getData.R`. 

Outputs:
- `luad_rna.csv`: a dataframe 
containing the barcodes of the RNAseq samples, the names of the corresponding files and
the corresponding patient IDs.
- `luad_mut_mc3.csv`: a dataframe containing mutation data for all patients  

### Data  filtering and pre-processing
1. To make list of patients containing all 3 data modalities (RNAseq, mutation and WSI), run `2_TCGA_getPatientList.R`
   
   Output:
   
   - `luad_patient_list.csv`: a dataframe with the IDs of the selected patients
   and their STK11 mutation status
     

2. To prepare the RNAseq counts for DE analysis, run `3_TCGA_prepareCountsTable.R`. This script combines all 
   RNAseq counts files into one single dataframe, for all patients in `luad_patient_list.csv`.
   Also allows the selection of gene types of interest (e.g. protein-coding genes). 
   
   Outputs:
   - `luad_rna_sample_info.csv`: dataframe containing RNAseq sample barcodes, file names,
   STK11 mutation status and sample type (tumor vs normal) for all patients in
     `luad_patient_list.csv`.
     
   - `luad_raw_counts_table.csv`: dataframe containing raw RNAseq counts 
     for all genes and all patients in `luad_patient_list.csv`
     
   - `luad_raw_counts_table_pc.csv`: dataframe containing raw RNAseq counts
   for protein-coding genes only, for all patients in `luad_patient_list.csv`

### DE analysis
To run the edgeR DE analysis pipeline, run `4_TCGA_runDEedgeR.R`.

Outputs:

- `sample_info_tumor_only_no_outliers.csv`: a filtered version of the
`luad_rna_sample_info.csv` dataframe, containing information for patients
  in `luad_patient_list.csv` only, and excluding normal tissue samples and
  outliers
  
- `luad_raw_counts_table_tumor_only_no_outliers.csv`: a filtered version of
`luad_raw_counts_table_pc.csv`, containing only the samples of interest
  
- `ExactTest_pc_genes_tumor_only_no_outliers.csv`: s dataframe with the DE
analysis results for all genes of interest (p-value, logFC, FDR...)
  
- `luad_rna_TMM_counts_pc_genes_tumor_only_no_outliers.csv`: TMM normalized
counts table for further analysis
  
- `genes.csv`: the list of gene symbols corresponding to the counts tables rows
  (order preserved)

### Data exploration
For dimensionality reduction and data visualization, run `5_TCGA_exploreData.R`

---

## Machine learning for feature selection
Method to derive a minimal gene expression signature for the classification of STK11
mutant tumors. After fitting a Logistic Regression or Random Forest model to the data, 
genes can be selected based on their feature importance.

![alt text](https://github.com/sophiehersz/TCGA_LUAD_project/blob/main/thumbnails/ML_pipeline.png?raw=true)

### Running the pipeline
Pipeline script: `nestedCV.py`

1. Data pre-processing - TMM normalized RNAseq counts are standardized (Z-score) and a min-max normalization is applied.
2. Creation and selection of gene lists of interest (feature sets) for model fitting
3. Fitting of a Logistic Regression and Random Forest model to the data using nested cross-validation
   1. inner loop - hyperparameter seach on training set (stratified 5-fold data split)
   2. outer loop - best model performance evaluation on test set (stratified 10-fold data split) using ROC-AUC metrics
4. Extraction of feature importance further for gene selection

Outputs:

- `TMM_counts_zscore.csv`: standardized counts table

- `TMM_counts_zscore_minmax.csv`: min-max normalized counts

---

## WSI data exploration

### Mean/median aggregation and dimensionality reduction
Pipeline script: `imageDataExploration.py`
1. Data pre-processing - computes mean or median embedding for each bag
2. Dimensionality reduction
3. Data visualization

Outputs:

- `mutation_status_embeddings.csv`: STK11 mutation status for WSI samples
- `mean_image_embeddings.csv`: mean embeddings for each WSI bag
- `median_image_embeddings.csv`: median embeddings for each WSI bag



### Classification using mean/median aggregation and Random Forest

Pipeline script: `imageDataClassification.py`

Uses mean or median aggregated image embeddings to train a Random Forest model
for STK11 mutation classification.

### Classification using Multiple-Instance Learning
Pipeline script: `imageDataMIL.py`

Attention-based Deep Multiple Instance Learning model (https://arxiv.org/abs/1802.04712)

1. Data filtering - selects files from patients in `sample_info_tumor_only_no_outliers.csv`
2. Data pre-processing
   1. Selects a random set of instances per bag (100 instances) to speed up training
   2. Selects a balanced subset of bags (same number of bags for each class)
3. Model training - trains an Attention-based Deep Multiple Instance Learning model
on training set, using stratified 5-fold cross-validation
4. Evaluates model on test set using ROC-AUC metrics

---

## Misc - Kaufman scoring
Pipeline script:`kaufmanScoring.py`

Based on publication *Kaufman et al., J Thorac Oncol, 2014*

1. Data pre-processing
   1. z-score standardization of TMM normalized counts
   2. min-max normalization
2. Data filtering - selection of genes of interest
3. Scoring - calculation of the unweighted mean of the standardized and
   normalized gene counts