library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)

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