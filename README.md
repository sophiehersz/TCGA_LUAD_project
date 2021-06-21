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
For RNAseq and mutation data download, run `1_TCGA_getData.R`

### Data  filtering and pre-processing
1. To make list of patients containing all 3 data modalities (RNAseq, mutation and WSI), run `2_TCGA_getPatientList.R`
2. To prepare the RNAseq counts for DE analysis, run `3_TCGA_prepareCountsTable.R`. This script combines all 
   RNAseq counts files into one single dataframe, for all patients belonging to the previous list.
   It also adds gene symbols to the dataframe and allows the selection of gene types of interest (e.g. protein-coding genes)

### DE analysis
To run the edgeR DE analysis pipeline, run `4_TCGA_runDEedgeR.R`.

### Data exploration
For dimensionality reduction and data visualization, run `5_TCGA_exploreData.R`

---

## Machine learning for feature selection
Method to derive a minimal gene expression signature for the classification of STK11
mutant tumors. After fitting a Logistic Regression or Random Forest model to the data, 
genes can be selected based on their feature importance.

![alt text](https://github.com/sophiehersz/TCGA_LUAD_project/blob/main/thumbnails/ML_pipeline.png?raw=true)

## Running the pipeline
Pipeline script: `nestedCV.py`

1. Data pre-processing - TMM normalized RNAseq counts are standardized (Z-score) and a min-max normalization is applied.
2. Selection of genes of interest (feature sets) for model fitting
3. Fitting of a Logistic Regression and Random Forest model to the data using nested cross-validation
   1. inner loop - hyperparameter seach on training set (stratified 5-fold data split)
   2. outer loop - best model performance evaluation on test set (stratified 10-fold data split)
4. Extraction of feature importance for gene selection