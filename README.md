# R scripts and processed datasets for classification of commercial cigarettes using inorganic, organic/isotopic, and low-level data fusion models

This repository contains the R scripts and processed datasets used to reproduce the chemometric analyses reported in the manuscript submitted to *Microchemical Journal*.

The repository includes three independent R workflows:

1. inorganic LDA model;
2. organic/isotopic LDA model;
3. low-level data fusion LDA model.

The scripts were prepared as supplementary reference material to ensure transparency, reproducibility, and traceability of the computational procedures used in the manuscript.

## Description

The study evaluates the classification of commercial cigarette samples using different analytical descriptor blocks and chemometric strategies.

The inorganic model uses elemental descriptors obtained by ICP-OES.

The organic/isotopic model uses n-alkane-related variables, bulk carbon isotope composition, carbon content, and nitrogen content.

The low-level data fusion model combines inorganic, organic, and isotopic descriptors at the variable level before chemometric modelling.

The main workflow includes:

- data preprocessing;
- log transformation, when applicable;
- autoscaling;
- block scaling for low-level data fusion;
- robust PCA-based multivariate outlier detection;
- removal of multivariate outliers;
- stratified train/test split;
- linear discriminant analysis;
- confusion matrix analysis;
- ROC/AUC analysis;
- variable-importance assessment;
- repeated train/test validation;
- permutation testing;
- robustness assessment.

## Repository structure

```text
low-level-data-fusion-cigarettes/
│
├── README.md
├── LICENSE
│
├── scripts/
│   ├── 01_inorganic_model_lda.R
│   ├── 02_organic_isotopic_model_lda.R
│   └── 03_low_level_data_fusion_robpca_lda.R
│
└── data/
    ├── inorganic_data.csv
    ├── pasta1.csv
    └── data_fusion_3.csv
```

## Scripts

### 01_inorganic_model_lda.R

This script reproduces the inorganic descriptor model.

It includes:

- selection of inorganic variables;
- descriptive statistics;
- ROBPCA-based multivariate outlier screening;
- multiclass LDA classification;
- binary LDA classification: Multinational vs Others;
- confusion matrices;
- ROC/AUC analysis;
- variable-importance analysis;
- repeated train/test validation;
- permutation testing.

The inorganic variables are:

```text
Zn, Sr, Ni, Mn, Fe, Cu, Co, Cd, Ba, B
```

### 02_organic_isotopic_model_lda.R

This script reproduces the organic/isotopic descriptor model.

It includes:

- selection of organic and isotopic variables;
- log transformation of n-alkane-related variables, when applicable;
- autoscaling;
- ROBPCA-based multivariate outlier screening;
- multiclass LDA classification;
- binary LDA classification: Multinational vs Others;
- confusion matrices;
- ROC/AUC analysis;
- variable-importance analysis;
- repeated train/test validation;
- permutation testing.

The organic/isotopic variables are:

```text
n-C18, n-C19, n-C21, n-C24, n-C29, n-C30, n-C31, n-C33,
C29/C31, ACL, alcanos totais, C18/total, Delta C13, %N, %C
```

### 03_low_level_data_fusion_robpca_lda.R

This script reproduces the low-level data fusion model.

It includes:

- selection of inorganic, organic, and isotopic descriptors;
- preprocessing of individual descriptor blocks;
- autoscaling by block;
- block scaling;
- low-level data fusion by variable concatenation;
- ROBPCA-based multivariate outlier screening;
- LDA classification;
- confusion matrices;
- ROC/AUC analysis;
- variable-importance analysis;
- comparison with individual descriptor blocks;
- repeated train/test validation;
- permutation testing;
- robustness assessment.

## Datasets

The folder `data/` contains the processed datasets used by the scripts.

### inorganic_data.csv

Processed dataset used for the inorganic LDA model.

### pasta1.csv

Processed dataset used for the organic/isotopic LDA model.

### data_fusion_3.csv

Processed dataset used for the low-level data fusion LDA model.

If a dataset is already loaded in the R environment with the same object name used in the original workflow, the scripts can also use the loaded object directly.

## How to reproduce the analyses

Open R or RStudio in the main repository folder.

To run the inorganic model:

```r
source("scripts/01_inorganic_model_lda.R")
```

To run the organic/isotopic model:

```r
source("scripts/02_organic_isotopic_model_lda.R")
```

To run the low-level data fusion model:

```r
source("scripts/03_low_level_data_fusion_robpca_lda.R")
```

Each script generates output files in the `outputs/` folder, including tables, figures, model summaries, confusion matrices, ROC/AUC results, variable-importance tables, robustness summaries, and session information.

## Required R packages

The scripts use the following R packages:

```text
tidyverse
caret
MASS
pROC
rrcov
ggplot2
ggrepel
```

Additional packages may be required depending on the specific script version and local R configuration.

If any required package is not installed, the scripts attempt to install it automatically.

## Output files

The scripts generate outputs such as:

- descriptive summary tables;
- ROBPCA outlier tables;
- training and test confusion matrices;
- LDA score tables;
- LDA score plots;
- explained variance tables;
- explained variance plots;
- ROC/AUC tables;
- ROC curves;
- variable-importance tables;
- variable-importance plots;
- repeated train/test validation summaries;
- permutation test summaries;
- session information files.

## Reproducibility

A fixed random seed is used in the scripts to improve reproducibility of train/test partitioning and robustness analyses.

The scripts also export session information files to document the R version and package versions used in the analyses.

## DOI

The archived version of this repository is available at Zenodo:

https://doi.org/10.5281/zenodo.20219595

DOI: 10.5281/zenodo.20219595

## Citation

If you use this repository, please cite the archived version available at Zenodo.

Suggested citation:

Rodrigues, L. S.; Massone, C. G.; Godoy, J. M. O. R scripts and processed datasets for classification of commercial cigarettes using inorganic, organic/isotopic, and low-level data fusion models. Version 1.0.0. Zenodo, 2026. DOI: 10.5281/zenodo.20219595.

## License

The R scripts are made available under the MIT License.

## Authors

Lucas Soares Rodrigues  
Carlos German Massone  
José Marcus de Oliveira Godoy











