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
```

## Files

The folder `scripts/` contains the R script used for the chemometric workflow.

The folder `data/` contains the processed dataset used in the analyses.

## How to reproduce the analyses

To reproduce the analyses, open R or RStudio and run:

```r
source("scripts/low_level_data_fusion_robpca_lda_publication.R")
```

The script will generate output files including confusion matrices, LDA scores, variable-importance tables, ROC/AUC results, robustness summaries, and figures.

## Required R packages

The script uses the following R packages:

- tidyverse
- caret
- MASS
- pROC
- rrcov
- ggplot2
- ggrepel

If any package is not installed, the script will attempt to install it automatically.

## Dataset

The processed dataset contains inorganic, organic, and isotopic descriptors obtained from commercial cigarette samples.

The inorganic descriptors include elemental concentrations determined by ICP-OES.

The organic and isotopic descriptors include n-alkane-related variables, bulk carbon isotope composition, carbon content, and nitrogen content.

## DOI

The archived version of this repository is available at Zenodo:

https://doi.org/10.5281/zenodo.20219595

DOI: 10.5281/zenodo.20219595

## Citation

If you use this repository, please cite the archived version available at Zenodo.

Suggested citation:

Rodrigues, L. S.; Massone, C. G.; Godoy, J. M. O. R scripts and processed dataset for low-level data fusion classification of commercial cigarettes. Version 1.0.0. Zenodo, 2026. DOI: 10.5281/zenodo.20219595.

## License

The R script is made available under the MIT License.

## Authors

Lucas Soares Rodrigues  
Carlos German Massone  
José Marcus de Oliveira Godoy
