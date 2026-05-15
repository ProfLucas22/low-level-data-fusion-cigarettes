# R scripts and processed dataset for low-level data fusion classification of commercial cigarettes

This repository contains the R script and processed dataset used to reproduce the chemometric analyses reported in the manuscript submitted to *Microchemical Journal*.

## Description

The workflow was developed for the classification of commercial cigarette samples using low-level data fusion of inorganic, organic, and isotopic descriptors.

The script performs the following steps:

- data preprocessing;
- low-level data fusion;
- block autoscaling;
- block scaling;
- robust PCA-based multivariate outlier detection;
- removal of multivariate outliers;
- stratified train/test split;
- linear discriminant analysis;
- confusion matrix analysis;
- ROC/AUC analysis;
- variable-importance assessment;
- repeated train/test validation;
- permutation testing;
- comparison of models with and without block scaling;
- evaluation of variable-importance stability.

## Repository structure

```text
low-level-data-fusion-cigarettes/
│
├── README.md
├── LICENSE
│
├── scripts/
│   └── low_level_data_fusion_robpca_lda_publication.R
│
└── data/
    └── data_fusion_3.csv