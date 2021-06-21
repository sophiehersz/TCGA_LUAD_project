# TCGA_LUAD_project

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

## RNAseq analysis pipeline
Differential expression (DE) analysis performed using edgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html)

![alt text](https://github.com/sophiehersz/TCGA_LUAD_project/blob/main/thumbnails/RNAseq_pipeline.png?raw=true)






