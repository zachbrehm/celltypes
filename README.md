Deconvolution Workflow: CIBERSORT
================
Zachary Brehm
3/3/2021

## Introduction

This repository outlines the steps of a workflow to perform
deconvolution of tissue gene expression data with
[CIBERSORT](https://cibersort.stanford.edu/). It requires reference
data, which we obtained with
[recount3](https://bioconductor.org/packages/release/bioc/html/recount3.html),
as well as mixture data. For the mixture data, we use coronary artery
samples from GTEx, also obtained with recount3.

## Preprocessing

This folder contains scripts to batch correct the reference data and
reduce the dimensions of these data for plotting. The batch correction
step is performed with
[RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html),
and the dimension reduction is shown with both the multidimensional
scaling and TSNE methods from with `cmdscale` and `Rtsne`. The following
plots serve as example results using each method.

![](images/gg_cmd.svg) ![](images/gg_tsne.svg)

## CIBERSORT

The cibersort directory contains the script ciber\_heat.R, used to
prepare the reference and mixture data for a deconvolution using the
CIBERSORT method linked above. There is also a script used to create a
heatmap displaying the proportion of the GTEx artery samples attributed
to each cell-type in the reference data. An example of this plot
follows.

![](images/gg_ciberHeat.png)
