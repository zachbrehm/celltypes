Deconvolution Workflow: CIBERSORT
================
Zachary Brehm
3/3/2021

## Introduction

This repository outlines the steps of a workflow to perform
deconvolution of tissue gene expression data with
[CIBERSORT](https://cibersort.stanford.edu/). Deconvolution methods are
used to gain insight into the composition of whole tissue samples. A
tissue is made of many different cell-types, and each cell-type
expresses unique genes indicative of their functions. Knowledge of how
much of each cell is present in a tissue can assist with differential
expression analyses.

CIBERSORT is one such method that estimates the proportion of each
cell-type found in a tissue. CIBERSORT requires 3 inputs.

-   **Reference** data, **R**<sub>*i* × *j*</sub>, a matrix of *i* genes
    by *j* samples of cell-type specific gene expression. Each cell-type
    that we are interested in estimating in the tissue should be present
    among these samples.

-   **Mixture** data, **M**<sub>*i* × *k*</sub>, a matrix of *i* genes
    by *k* samples of whole tissue gene expression. These are the
    samples we will deconvolute.

-   A **phenotype** indication table, **P**<sub>*h* × *j*</sub>, a
    matrix of *h* cell-types by *j* samples. The columns here must match
    to the samples in the reference data. Here, entry
    *p*<sub>*m**n*</sub> = 1 if sample *n* is of cell-type *m*,
    otherwise it will equal 2.

All data used for this example was obtained with
[recount3](https://bioconductor.org/packages/release/bioc/html/recount3.html).
For the mixture data, we use coronary artery samples from
[GTEx](https://www.gtexportal.org/home/), and the reference samples are
from cell-types found in atherosclerotic artery. For our purposes, we
will use the following cell-types:

-   Smooth muscle cells
-   Endothelial cells
-   Erythrocytes
-   Macrophages
-   Lymphocytes
-   Adipocytes
-   Cardiac muscle cells

## Preprocessing

This folder contains scripts to batch correct the reference data and
reduce the dimensions of these data for plotting. Since our reference
data are from numerous different experiments, The batch correction step
is performed with `RUVr` from the
[RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
package, and the dimension reduction is shown with both the
multidimensional scaling and TSNE methods from with `cmdscale` in base R
and `Rtsne` from the `Rtsne` package. The following plots serve as
example results using each method.

![](images/gg_cmd.svg) ![](images/gg_tsne.svg)

## CIBERSORT

The cibersort directory contains the script ciber\_heat.R, used to
prepare the reference and mixture data for a deconvolution using the
CIBERSORT method linked above. There is also a script used to create a
heatmap displaying the proportion of the GTEx artery samples attributed
to each cell-type in the reference data. An example of this plot
follows.

![](images/gg_ciberHeat.png)
