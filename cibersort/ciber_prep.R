library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)

## cibersort takes in a reference expression set to deconvolute a mixture expression set
## read in the reference set
rse <- readRDS(file = "data/ruv_reference_2021_v2.Rds")

## read in the mixture set
gtex <- readRDS(file = "data/rse_coronary_recount3.Rds")

## deseq requires the first assay to be counts
## recount3 assigns raw counts to the first slot, so we
## remove raw counts assay from mixture data leaving only the transformed counts
assays(gtex)[1] <- NULL

## filter low count genes to ease computation
lowCounts <- apply(assay(gtex), 1, function(x){all(x < 1000)})
gtex <- gtex[-which(lowCounts), ]

## identify the indices of celltypes of each sample in the reference set
celltypes <- unique(rse$celltype)
celltypeInd <- list(cmc = vector(mode = "numeric", length = sum(cell)))
cardiacInd <- which(rse$celltype == "CardiacMuscleCell")
endoInd <- which(rse$celltype == "EndothelialCell")
lymphInd <- which(rse$celltype == "Lymphocytes")
macroInd <- which(rse$celltype == "Macrophage")
erythInd <- which(rse$celltype == "Erythrocyte")
fatInd <- which(rse$celltype == "FatCell")
smcInd <- which(rse$celltype == "SmoothMuscleCell")

## cibersort requires a phenotype matrix as an input, nCelltypes by nSamples
## the required format uses the value 1 in entry i j to indicate that
## sample j is of celltype i, otherwise the entry takes on the value of 2
## create phenotype identification matrix 
classMat <- matrix(data = rep.int(x = 2, times = length(celltypes) * ncol(rse)), 
                   nrow = length(celltypes), ncol = ncol(rse))
colnames(classMat) <- colnames(rse)

## mark positions of each celltype in phenotype matrix with 1
classMat[1, cardiacInd] <- 1
classMat[2, endoInd] <- 1
classMat[3, lymphInd] <- 1
classMat[4, smcInd] <- 1
classMat[5, macroInd] <- 1
classMat[6, erythInd] <- 1
classMat[7, fatInd] <- 1

## reformat the phenotype matrix according to cibersort specifications
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
write.table(x = classes, file = "cibersort/input/ciber_phe_recount3.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(x = reference, file = "cibersort/input/ciber_ref_recount3.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(x = mixture, file = "cibersort/input/ciber_mix_recount3.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = TRUE)