library(DESeq2)
library(tidyverse)
library(viridis)

rse <- readRDS("ruv_rse_noFoamFibro_k20.Rds")
rse <- rse[-which(low_counts), ]
dds <- DESeqDataSet(rse, design = ~1)
vst <- varianceStabilizingTransformation(dds)

celltypes <- levels(vst$celltype)

cellInds <- vector(mode = "list", length = length(celltypes))
names(cellInds) <- celltypes

dfList <- vector(mode = "list", length = length(celltypes))
names(dfList) <- celltypes

for(i in 1:length(celltypes)){
  cellInds[[i]] <- which(vst$celltype == celltypes[i])
  dfList[[i]] <- data.frame(Gene = rownames(vst), 
                            Mean = rowMeans(assay(vst)[,cellInds[[i]]]),
                            SD = rowVars(assay(vst)[,cellInds[[i]]]) %>% sqrt())
}


colorPal <- viridis(3, begin = 0.25, end = 0.75)
plots <- lapply(dfList, function(x){ggplot(x, aes(x = Mean, y = SD)) + geom_point(alpha = 0.5) +
    geom_hline(yintercept = quantile(x$SD[x$Mean >= 5], 0.9), color = colorPal[1], show.legend = TRUE) + 
    geom_text(aes(x = 0, y = quantile(x$SD[x$Mean>= 5], 0.9), label = paste("90:", round(quantile(x$SD[x$Mean>= 5], 0.9), digits = 3)), vjust = -0.1)) +
    geom_hline(yintercept = quantile(x$SD[x$Mean >= 5], 0.95), color = colorPal[2]) + 
    geom_text(aes(x = 1.25, y = quantile(x$SD[x$Mean >= 5], 0.95), label = paste("95:", round(quantile(x$SD[x$Mean >= 5], 0.95), digits = 3)), vjust = -0.1)) +
    geom_hline(yintercept = quantile(x$SD[x$Mean >= 5], 0.99), color = colorPal[3]) +
    geom_text(aes(x = 2.5, y = quantile(x$SD[x$Mean >= 5], 0.99), label = paste("99:", round(quantile(x$SD[x$Mean >= 5], 0.99), digits = 3)), vjust = -0.1)) +
    geom_hline(yintercept = quantile(x$SD[x$Mean >= 5], 0.1), color = colorPal[1]) + 
    geom_text(aes(x = 0, y = quantile(x$SD[x$Mean >= 5], 0.1), label = paste("10:", round(quantile(x$SD[x$Mean >= 5], 0.10), digits = 3)), vjust = -0.1)) +
    geom_hline(yintercept = quantile(x$SD[x$Mean >= 5], 0.05), color = colorPal[2]) + 
    geom_text(aes(x = 1.25, y = quantile(x$SD[x$Mean >= 5], 0.05), label = paste("5:", round(quantile(x$SD[x$Mean >= 5], 0.05), digits = 3)), vjust = -0.1)) +
    geom_hline(yintercept = quantile(x$SD[x$Mean >= 5], 0.01), color = colorPal[3]) +
    geom_text(aes(x = 2.5, y = quantile(x$SD[x$Mean >= 5], 0.01), label = paste("1:", round(quantile(x$SD[x$Mean >= 5], 0.01), digits = 3)), vjust = -0.1))})
for(i in 1:length(celltypes)){
  plots[[i]] <- plots[[i]] + labs(title = paste("Mean vs SD for VST expression per gene from", celltypes[i]))
}
varHist <- vector(mode = "list", length = length(celltypes))

varHist[[1]] <- ggplot(varDf, aes(x = CardiacMuscleCell)) + geom_histogram(bins = 50)
varHist[[2]] <- ggplot(varDf, aes(x = EndothelialCell)) + geom_histogram(bins = 50)
varHist[[3]] <- ggplot(varDf, aes(x = Erythrocyte)) + geom_histogram(bins = 50)
varHist[[4]] <- ggplot(varDf, aes(x = FatCell)) + geom_histogram(bins = 50)
varHist[[5]] <- ggplot(varDf, aes(x = Fibroblast)) + geom_histogram(bins = 50)
varHist[[6]] <- ggplot(varDf, aes(x = Foam)) + geom_histogram(bins = 50)
varHist[[7]] <- ggplot(varDf, aes(x = Lymphocytes)) + geom_histogram(bins = 50)
varHist[[8]] <- ggplot(varDf, aes(x = Macrophage)) + geom_histogram(bins = 50)
varHist[[9]] <- ggplot(varDf, aes(x = SmoothMuscleCell)) + geom_histogram(bins = 50)

erythGenes_lowSd <- which(dfList[[3]]$Mean >= 20 & dfList[[3]]$SD <= 0.23)
erythHighMean_lowSD_gtex <- t(assay(vst_coronary)[erythGenes_lowSd,])
colnames(erythHighMean_lowSD_gtex) <- c("GYPC", "RN7SL4P", "RN7SL2", "RN7SL1", "Not.In.ENSMBL")
erythHighMean_lowSD_gtex <- erythHighMean_lowSD_gtex[,sort(dfList[[3]]$SD[erythGenes_lowSd], decreasing = TRUE, index.return = TRUE)$ix]
erythHighMean_lowSD_rse <- apply(assay(vst)[erythGenes_lowSd,vst$celltype == "Erythrocyte"], 1, summary)
colnames(erythHighMean_lowSD_rse) <- c("GYPC", "RN7SL4P", "RN7SL2", "RN7SL1", "Not.In.ENSMBL")
erythHighMean_lowSD_rse <- erythHighMean_lowSD_rse[,sort(dfList[[3]]$SD[erythGenes_lowSd], decreasing = TRUE, index.return = TRUE)$ix]

erythGenes_highSd <- which(dfList[[3]]$Mean >= 15 & dfList[[3]]$SD >= 1.421)
erythHighMean_highSD_gtex <- apply(assay(vst_coronary)[erythGenes_highSd,], 1, summary)
colnames(erythHighMean_highSD_gtex) <- c("SLC4A1", "E2F2", "MAP2K3", "ATP5F1E", "TERF2", "MBOAT2", "ZEB1", "KRT1", "IKZF1", "HBA2", "RNY1",
                                         "RPL12P16", "HBA1", "HBB", " AC092490.1")
erythHighMean_highSD_gtex <- erythHighMean_highSD_gtex[,sort(dfList[[3]]$SD[erythGenes_highSd], decreasing = TRUE, index.return = TRUE)$ix]
erythHighMean_highSD_rse <- apply(assay(vst)[erythGenes_highSd,vst$celltype == "Erythrocyte"], 1, summary)
colnames(erythHighMean_highSD_rse) <- hgncVec[erythGenes_highSd]
erythHighMean_highSD_rse <- erythHighMean_highSD_rse[,sort(dfList[[3]]$SD[erythGenes_highSd], decreasing = TRUE, index.return = TRUE)$ix]

write.csv(erythHighMean_highSD_gtex, file = "erythHighMean_highSD_gtex.csv")
write.csv(erythHighMean_highSD_rse, file = "erythHighMean_highSD_rse.csv")
write.csv(erythHighMean_lowSD_gtex, file = "erythHighMean_lowSD_gtex.csv")
write.csv(erythHighMean_lowSD_rse, file = "erythHighMean_lowSD_rse.csv")

eryth1 <- sort(assay(vst)[,which(vst$run == "SRR2124301")], decreasing = TRUE, index.return = TRUE)
eryth2 <- sort(assay(vst)[,which(vst$run == "SRR2124300")], decreasing = TRUE, index.return = TRUE)
eryth3 <- sort(assay(vst)[,which(vst$run == "SRR2124299")], decreasing = TRUE, index.return = TRUE)
eryth4 <- sort(assay(vst)[,which(vst$run == "SRR2038798")], decreasing = TRUE, index.return = TRUE)
top50GenesErythrocytes <- data.frame(SRR2124301 = hgnc[eryth1$ix[1:50]], 
                                     SRR2124300 = hgnc[eryth2$ix[1:50]],
                                     SRR2124299 = hgnc[eryth3$ix[1:50]],
                                     SRR2038798 = hgnc[eryth4$ix[1:50]])
