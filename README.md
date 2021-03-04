Deconvolution Workflow: CIBERSORT
================
Zachary Brehm
2021-03-04

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
    *p*<sub>*m*, *n*</sub> = 1 if sample *n* is of cell-type *m*,
    otherwise it will equal 2.

All data used for this example was obtained with
[recount3](https://bioconductor.org/packages/release/bioc/html/recount3.html).
For the mixture data, we use coronary artery samples from
[GTEx](https://www.gtexportal.org/home/), and the reference samples are
from cell-types found in atherosclerotic artery. Atherosclerosis is a
disease which is strongly characterized by compositional changes in
tissue. As this disease progresses, plaques develop along the inside of
the artery, causing health complications or death. The nature of this
disease makes it a natural candidate for deconvolution. For example, in
a differential expression analyses, we would like to investigate whether
shifts in gene expression can be attributed to differences in form
between samples or differences in function between samples. For our
purposes, we will use the following cell-types:

-   Smooth muscle cells
-   Endothelial cells
-   Erythrocytes
-   Macrophages
-   Lymphocytes
-   Adipocytes
-   Cardiac muscle cells

## Downloading the data

To use `recount3`, you must first install Bioconductor. This can be
accomplished by running

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
```

We will also need `DESeq2`, `RUVSeq`, and `edgeR` later on, so we can
install all of these at once. These must be installed with
`BiocManager`. After installing Bioconductor, run

``` r
BiocManager::install(c("DESeq2", "edgeR", "RUVSeq", "recount3"))
```

One source for reference data is the Sequence Read Archive
([SRA](https://www.ncbi.nlm.nih.gov/sra)), much of which is available
through `recount3` along with data from GTEx. For example, we can
download the Blood Vessel tissue data from GTEx and subset to only the
samples containing coronary artery.

``` r
## downloading gtex data with recount3
library(recount3)

## get table of available projects from recount3
human_projects <- available_projects()

## find gtex projects
gtex <- subset(human_projects, file_source == "gtex" & project_type == "data_sources")

## gtex data in recount3 is stored according to tissue type
## we want the blood vessel data in this case
proj_info <- subset(gtex, project == "BLOOD_VESSEL" & project_type == "data_sources")

## make the summarized experiment object from the blood vessel gtex data
rse_blood_vessel <- create_rse(proj_info)

## we only want the coronary artery samples, stored in gtex.smtsd as "Artery - Coronary"
rse_coronary <- rse_blood_vessel[ , rse_blood_vessel$gtex.smtsd == "Artery - Coronary"]

## transform counts as with the reference data
assay(rse_coronary, "counts") <- transform_counts(rse_coronary)
saveRDS(rse_coronary, file = "data/rse_coronary.Rds")
```

## Preprocessing

This folder contains scripts to batch correct the reference data and
reduce the dimensions of these data for plotting. Since our reference
data are from numerous different experiments, this step estimates and
removes technical variation so we can better represent the underlying
biology that we wish to estimate. The batch correction step is performed
with `RUVr` from the
[RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
package. This script takes in the reference data that was acquired with
`recount3` and outputs a batch corrected version of the summarized
experiment object.

``` r
## read in reference data for batch correction
rse <- readRDS(file = "data/rse_reference.Rds")

## filter out low counts to ease computation
low_counts <- apply(assay(rse, "counts"), MARGIN = 1, FUN = function(x){all(x < 1000)})
rse <- rse[-which(low_counts), ]

## create expression set for edgeR
celltypes <- rse$celltype
set <- newSeqExpressionSet(assay(rse, "counts"), 
                           phenoData = data.frame(celltypes, row.names = colnames(rse)))
design <- model.matrix(~celltypes, data = pData(set))
genes <- rownames(rse)

## run edgeR
y <- DGEList(counts = counts(set), group = celltypes) %>% 
  calcNormFactors(object = ., method = "TMM") %>%
  estimateGLMCommonDisp(y = ., design = design) %>%
  estimateGLMTagwiseDisp(y = ., design = design)
fit <- glmFit(y = y, design = design)
res <- residuals(fit, type = "pearson")

## run RUV seq with RUVr method for blind batch correction
set0 <- RUVr(set, genes, k = 20, res)

## create summarized experiment object for cibersort preparation
## see ciber_prep.R for further steps
ruv <- SummarizedExperiment(assays = list(counts = set0@assayData$normalizedCounts), 
                            rowData = rowData(rse), colData = colData(rse))
saveRDS(ruv, file = "data/ruv_reference.Rds")
```

We also want to plot the data to see if the samples from each cell-type
are reasonably similar to one another. To do so, we calculate the
pairwise euclidean distance between samples. This distance matrix is
then used for a dimension reduction step to allow for plotting in a
plane, and is demonstrated with both the multidimensional scaling and
TSNE methods from with `cmdscale` in base R and `Rtsne` from the `Rtsne`
package. The following plots serve as example results using each method.
These are presented as static images here, but with the `plotly`
package, we can create interactive plots and mouse over the points to
identify which study the samples come from. To do so, simply store the
`ggplot2` output and feed it into the `ggplotly` function.

``` r
## read in batch corrected gene expression data
rse <- readRDS(file = "data/ruv_reference.Rds")

## make rse into dds object for variance stabilizing
dds <- DESeqDataSet(rse, design = ~1)

## create vst object to attenuate the mean/variance relationship
## of the expression data
vst <- varianceStabilizingTransformation(dds)

## compute euclidean distance between samples
dist <- dist(t(assay(vst)))

## compute multidimensional scaling of distance matrix
## and convert to dataframe to plot
cmdDf <- cmdscale(dist) %>% 
  data.frame(.) %>% 
  rownames_to_column(var = "Sample") %>% 
  mutate(Celltype = rse$celltype, 
         Project = rse$study) %>% 
  select(Sample, 
         Celltype,
         Project,
         everything())

gg_cmd <- ggplot(cmDf, aes(x = V1, y = V2, color = Celltype, alpha = 0.8, text = Project)) + 
  geom_point() + 
  labs(title = "Dimension Reduction with CMD", 
       x = "First Dimension", 
       y = "Second Dimension", 
       color = "Cell-Type") + 
  scale_color_viridis_d(labels = c("Cardiac Muscle", 
                                   "Endothelial", 
                                   "Erythrocyte", 
                                   "Adipocyte", 
                                   "Lymphocyte", 
                                   "Macrophage", 
                                   "Smooth Muscle"),
                        begin = 0, end = 0.95) + 
  guides(alpha = FALSE) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line("grey30", linetype = "dotted"),
        panel.grid.minor = element_line("grey30", linetype = "dotted"),
        legend.key = element_rect(fill = "white"))

ggsave(filename = "preprocessing/gg_cmd.svg", plot = gg_cmd, 
       device = "svg", dpi = 600, width = 11, height = 8.5, units = "in")
```

![](images/gg_cmd.svg)

With `Rtsne`:

``` r
## run tsne with the previously created distance matrix
cellTsne <- Rtsne(dist, is_distance = TRUE, perplexity = 40)

## make tsne dataframe and plot
tsneDf <- data.frame(X = cellTsne$Y[,1], 
                     Y = cellTsne$Y[,2], 
                     Z = rse$celltype,
                     Project = rse$study)

gg_tsne <- ggplot(tsneDf, aes(x = X, y = Y, color = Z, text = Project, alpha = 0.8)) + 
  geom_point() + 
  labs(title = "Dimension Reduction with TSNE", 
       x = "First Dimension", 
       y = "Second Dimension", 
       color = "Cell-Type",
       text = rse$study) + 
  guides(alpha = FALSE) + 
  scale_color_viridis_d(labels = c("Cardiac Muscle", 
                                   "Endothelial", 
                                   "Erythrocyte", 
                                   "Adipocyte", 
                                   "Lymphocyte", 
                                   "Macrophage", 
                                   "Smooth Muscle"),
                        begin = 0, end = 0.95) + 
                        ## with end = 1 the yellow is a bit too bright for a white background
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line("grey30", linetype = "dotted"),
        panel.grid.minor = element_line("grey30", linetype = "dotted"),
        legend.key = element_rect(fill = "white"))

ggsave(filename = "preprocessing/gg_tsne.svg", plot = gg_tsne, 
       device = "svg", dpi = 600, width = 10.5, height = 8, units = "in")
```

![](images/gg_tsne.svg)

## CIBERSORT

The cibersort directory contains the script ciber\_heat.R, used to
prepare the reference and mixture data for a deconvolution using the
CIBERSORT method linked above.

``` r
## cibersort takes in a reference expression set to deconvolute a mixture expression set
## read in the reference set
rse <- readRDS(file = "data/ruv_reference.Rds")

## read in the mixture set
gtex <- readRDS(file = "data/rse_coronary.Rds")

## deseq requires the first assay to be counts
## recount3 assigns raw counts to the first slot, so we
## remove raw counts assay from mixture data leaving only the transformed counts
assays(gtex)[1] <- NULL

## filter low count genes to ease computation
lowCounts <- apply(assay(gtex), 1, function(x){all(x < 1000)})
gtex <- gtex[-which(lowCounts), ]

## identify the indices of celltypes of each sample in the reference set
celltypes <- unique(rse$celltype)
celltypeInd <- lapply(celltypes, FUN = function(x) {which(rse$celltype == x)})
names(celltypeInd) <- celltypes

## cibersort requires a phenotype matrix as an input, nCelltypes by nSamples
## the required format uses the value 1 in entry i j to indicate that
## sample j is of celltype i, otherwise the entry takes on the value of 2
## create phenotype identification matrix 
classMat <- matrix(data = rep.int(x = 2, times = length(celltypes) * ncol(rse)), 
                   nrow = length(celltypes), ncol = ncol(rse))
colnames(classMat) <- colnames(rse)

## mark positions of each celltype in phenotype matrix with 1
for(i in 1:length(celltypes)){
  classMat[i,celltypeInd[[i]]] <- 1
}

## reformat the phenotype matrix according to cibersort specifications
## abbreviating names of cells for easier formatting
classes <- as.data.frame(classMat) %>%
  mutate(Class = c("CMC", "END", "LYM", "SMC", "MAC", "RBC", "ADI")) %>%
  select(Class, everything())

## normalize the reference set
ref_dds <- DESeqDataSet(rse, design = ~1) %>% estimateSizeFactors()

## normalize the mixture set
mix_dds <- DESeqDataSet(gtex, design = ~1) %>% estimateSizeFactors()

## transform reference data into the format required by cibersort
## first column gene symbols, remaining columns are normalized count data
reference <- counts(ref_dds, normalized = TRUE) %>% round() %>% 
  as.data.frame() %>% mutate(HGNC = rowData(ref_dds)$gene_name) %>% select(HGNC, everything())

## transform mixture data into the format required by cibersort
## first column gene symbols, remaining columns are normalized count data
mixture <- counts(mix_dds, normalized = TRUE) %>% round() %>% 
  as.data.frame() %>% mutate(HGNC = rowData(mix_dds)$gene_name) %>% select(HGNC, everything())

## write final tables to input into cibersort
write.table(x = classes, file = "cibersort/input/ciber_phe.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(x = reference, file = "cibersort/input/ciber_ref.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(x = mixture, file = "cibersort/input/ciber_mix.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = TRUE)
```

There is also a script used to create a heatmap displaying the
proportion of the GTEx artery samples attributed to each cell-type in
the reference data. This can be made from the output that you would
receive from CIBERSORT.

``` r
## read in proportion estimates from cibersort
ciberOut <- read.table(file = "cibersort/output/ciberOut_gtexCoronary.txt", header = TRUE)

## extract sample and proportion values and pivot to long format for plotting
res_ciber_df <- data.frame(ciberOut[,1:8]) %>% 
  pivot_longer(cols = 2:8, names_to = "Celltype", values_to = "Proportion")

## plot proportion estimates in a heat map, rows indexed by samples and
## columns indexed by cell types
## color filled with estimated proportion of tissue of composition per cell
ggplot(res_ciber_df, 
       aes(x = Celltype, y = InputSample, fill = Proportion)) + geom_tile() + 
  scale_fill_viridis_c(breaks = seq(0, 1, 0.1), 
                       limits = c(0, 1), 
                       label = seq(0, 100, 10)) + 
  labs(title = "Reference Based Deconvolution of GTEx Coronary Artery with CIBERSORT",
       x = "Cell-Type",
       y = "Sample",
       fill = "% of Tissue") + 
  theme(text = element_text(family = "Fira Sans"),
        plot.background = element_rect(fill = "white"), 
        axis.text.x = element_text(angle = 45, vjust = 0.65), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key.size = unit(3,"line"))
```

![](images/gg_ciberHeat.png)

-   Smooth muscle cells are the most abundant cell-type. Smooth muscle
    cells are the primary cell-type found in coronary artery, as they
    are a main structural component of the artery itself.

-   Endothelial cells line the interior of the artery, so they are
    consistently found in small amounts across all samples.

-   Erythrocytes are essentially nonexistant, which we expect to see
    given that these samples are washed so little blood should remain in
    them.

-   Macrophages and lymphocytes are the main indicators of inflammation
    and disease in these samples, weighted towards macrophages in
    particular. They appear in small amounts more less consistently than
    endothelial cells.

-   Adipose tissue and cardiac muscle are not apart of a typical
    coronary artery, but appear here since the GTEx samples present with
    some excess tissue attached that was not trimmed off when harvested.
    This is made apparent by the sporadic and sometimes large amounts of
    adipose signal in the heatmap. Cardiac muscle appears in small
    amounts more consistently, but as muscle tissue it shares some
    similar characteristics with smooth muscle, making it difficult to
    completely isolate these cell-types away from one another.
