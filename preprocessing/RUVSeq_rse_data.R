library(tidyverse)
library(edgeR)
library(RUVSeq)

rse <- readRDS(file = "data/rse_reference_2021_v1.Rds")

assays(rse)[1] <- NULL
#filterGenes <- which(rownames(rse) %in% rownames(vst))

low_counts <- apply(assay(rse, "counts"), MARGIN = 1, FUN = function(x){all(x < 1000)})
rse <- rse[-which(low_counts), ]



#zeros <- which(colSums(assay(rse,"counts")) == 0)
#rse <- rse[,-zeros]
celltypes <- rse$celltype
set <- newSeqExpressionSet(assay(rse, "counts"), 
                           phenoData = data.frame(celltypes, row.names = colnames(rse)))
design <- model.matrix(~celltypes, data = pData(set))
genes <- rownames(rse)
y <- DGEList(counts = counts(set), group = celltypes) %>% 
  calcNormFactors(object = ., method = "TMM") %>%
  estimateGLMCommonDisp(y = ., design = design) %>%
  estimateGLMTagwiseDisp(y = ., design = design)

fit <- glmFit(y = y, design = design)
res <- residuals(fit, type = "pearson")

set0 <- RUVr(set, genes, k = 20, res)

ruv <- SummarizedExperiment(assays = list(counts = set0@assayData$normalizedCounts), 
                            rowData = rowData(rse), colData = colData(rse))
saveRDS(ruv, file = "data/ruv_reference_2021_v1.Rds")

#rm(fit, rse, ruv, y, res)
tmprse <- SummarizedExperiment(assays = assay(rse, "counts"), colData = colData(rse), rowRanges = rowRanges(rse))
dds <- DESeq2::DESeqDataSet(tmprse[,-which(colSums(assay(tmprse)) == 0)], design = ~1) %>% DESeq2::estimateSizeFactors() %>% DESeq2::estimateDispersions()
