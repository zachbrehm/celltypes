library(tidyverse)
library(edgeR)
library(RUVSeq)

## read in reference data for batch correction
rse <- readRDS(file = "data/rse_reference_2021_v1.Rds")

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
saveRDS(ruv, file = "data/ruv_reference_2021_v2.Rds")
