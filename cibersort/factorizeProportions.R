library(DESeq2)
library(tidyverse)
library(recount)
ciberOut <- ciberOut <- read.table(file = "cibersort/output/ciberOut_recount3.txt", header = TRUE)
celltypeProportions <- ciberOut[,1:8]

gtex <- readRDS(file = "data/rse_coronary_recount3.Rds")
assays(gtex)[1] <- NULL
cac <- read_excel("Coronary_artery_case_images_8.27.20.xlsx")

gtexNames <- str_split(colnames(gtex), "-") %>% sapply(FUN = function(x) paste(x[1], x[2], x[3], sep = "-"))

cacSub <- subset(cac, `Specimen ID` %in% gtexNames)
cacSub <- cacSub[-148,]

duplicatedSamples <- which(duplicated(gtexNames))

gtex <- gtex[,-c(49, 88, 212, 213, 214)]
gtexNames <- gtexNames[-c(49, 88, 212, 213, 214)]

map <- match(gtexNames, cacSub$`Specimen ID`)

cacSub <- cacSub[map,]

lowCounts <- apply(assay(gtex), 1, function(x){all(x < 1000)})
gtex <- gtex[-which(lowCounts), ]

proportionFactors <- celltypeProportions
celltypes <- colnames(celltypeProportions)[2:8]
proportionFactors$CMC <- factor((celltypeProportions$CMC >= quantile(celltypeProportions$CMC, seq(0, 1, 0.05))['95%'])*1)
proportionFactors$LYM <- factor((celltypeProportions$LYM > 0)*1)
proportionFactors$ADI[km$cluster == 2] <- 0
proportionFactors$ADI[km$cluster == 1] <- 1
proportionFactors$ADI[km$cluster == 3] <- 2
proportionFactors$ADI <- factor(proportionFactors$ADI)

proportionFactors$END <- round(proportionFactors$END, digits = 4)
proportionFactors$SMC <- round(proportionFactors$SMC, digits = 4)
proportionFactors$MAC <- round(proportionFactors$MAC, digits = 4)

proportionFactors <- proportionFactors[-c(49, 88, 212, 213, 214),]

gtex$Plaque <- factor(cacSub$`Plaque Type`)
gtex$PlaqueInd <- factor((cacSub$`Plaque Type` > 0)*1)
gtex$CMC <- proportionFactors$CMC
gtex$END <- proportionFactors$END
gtex$LYM <- proportionFactors$LYM
gtex$SMC <- proportionFactors$SMC
gtex$MAC <- proportionFactors$MAC
gtex$ADI <- proportionFactors$ADI

deseqPlaque <- DESeqDataSet(gtex, design = ~ Plaque) %>% DESeq()
deseqPlaqueInd <- DESeqDataSet(gtex, design = ~ PlaqueInd) %>% DESeq()
deseqPlaque_cells <- DESeqDataSet(gtex, design = ~ Plaque + CMC + END + LYM + SMC + MAC + ADI) %>% DESeq()
deseqPlaqueInd_cells <- DESeqDataSet(gtex, design = ~ PlaqueInd + CMC + END + LYM + SMC + MAC + ADI) %>% DESeq()

resPlaque <- results(deseqPlaque, contrast = c("Plaque", 4, 0))
resPlaque$HGNC <- rowData(gtex)$gene_name
resPlaque <- resPlaque[order(resPlaque$padj, decreasing = FALSE), c(7, 1:6)]
resPlaque_filter <- resPlaque[resPlaque$baseMean > 500, ]
resPlaque_filter$padj <- p.adjust(resPlaque_filter$pvalue, method = "BH")

resPlaqueInd <- results(deseqPlaqueInd, contrast = c("PlaqueInd", 1, 0))
resPlaqueInd$HGNC <- rowData(gtex)$gene_name
resPlaqueInd <- resPlaqueInd[order(resPlaqueInd$padj, decreasing = FALSE), c(7, 1:6)]
resPlaqueInd_filter <- resPlaqueInd[resPlaqueInd$baseMean > 500,]
resPlaqueInd_filter$padj <- p.adjust(resPlaqueInd_filter$pvalue, method = "BH")

resPlaque_cells <- results(deseqPlaque_cells, contrast = c("Plaque", 4, 0))
resPlaque_cells$HGNC <- rowData(gtex)$gene_name
resPlaque_cells <- resPlaque_cells[order(resPlaque_cells$padj, decreasing = FALSE), c(7, 1:6)]
resPlaque_cells_filter <- resPlaque_cells[resPlaque_cells$baseMean > 500, ]
resPlaque_cells_filter$padj <- p.adjust(resPlaque_cells_filter$pvalue, method = "BH")

resPlaqueInd_cells <- results(deseqPlaqueInd_cells, contrast = c("PlaqueInd", 1, 0))
resPlaqueInd_cells$HGNC <- rowData(gtex)$gene_name
resPlaqueInd_cells <- resPlaqueInd_cells[order(resPlaqueInd_cells$padj, decreasing = FALSE), c(7, 1:6)]
resPlaqueInd_cells_filter <- resPlaqueInd_cells[resPlaqueInd_cells$baseMean > 500, ]
resPlaqueInd_cells_filter$padj <- p.adjust(resPlaqueInd_cells_filter$pvalue, method = "BH")

celltypeGenes <- list(TCells = c("CCL5", "CD8A", "GZMA", "GZMK", 
                                "LTB", "MAL", "IL32", "IL7R",
                                "FGFBP2", "GZMB", "PRF1", "ADGRG1",
                                "CD2", "CD69", "SPOCK2", "CD3D"),
                     BCells = c("CD79A", "BANK1", "BCL11A", "CD19"),
                     Myeloid = c("C1QA", "C1QB", "C1QC", "MS4A6A",
                                 "S100A8", "S100A9", "SERPINA1", "FCN1",
                                 "MMP9", "SPP1", "GPNMB", "APOC1",
                                 "CLEC10A", "FCER1A", "HLA-DQA1", "CPVL"),
                     Mast = c("CMA1", "HDC", "KIT", "TPSAB1"),
                     Mixed = c("LOC100131257", "PGM5P2", "MAB21L3", "TMEM212"),
                     SMC = c("AEBP1", "TAGLN", "MYH11", "PDGFRB", 
                             "NOTCH3", "MFAP4", "ACTA2", "MGP",
                             "COL1A1", "COL1A2", "COL3A1"),
                     END = c("COL4A1", "COL4A2", "SPARCL1", "PLVAP",
                             "MPZL2", "SULF1", "VWF", "EDN1")
)

atheroGenes <- c("CDKN2A", "CDKN2B", "MTHFD1L", "CELSR2", "PSRC1", "SORT1", "MIA3", "SMAD3",
                 "APOE", "APOC1", "APOC4", "LDLR", "APOB", "PCKS9", "CHRNA3")

genes <- unlist(celltypeGenes, use.names = FALSE)
genes <- c(genes[-which(genes == "APOC1")], atheroGenes)
geneSub <- genes[which(genes %in% resPlaque$HGNC)]
genesDouble <- geneSub %>% rep(times = 4)

asso <- c(rep("TCell", times = 3),
          rep("Myeloid", times = 11),
          rep("Mast", times = 1),
          rep("SMC", times = 11),
          rep("END", times = 8),
          rep("Athero", times = 7)
          ) %>% rep(times = 4)
adjusted <- factor(rep(c("No", "Yes"), each = 65)) %>% rep(.,times = 2)

plaque <- factor(rep(c("Allv0", "4v0"), each = 130))

map1 <- match(geneSub, resPlaqueInd_filter$HGNC[resPlaqueInd_filter$HGNC %in% geneSub])
map2 <- match(geneSub, resPlaqueInd_cells_filter$HGNC[resPlaqueInd_cells_filter$HGNC %in% geneSub])
map3 <- match(geneSub, resPlaque_filter$HGNC[resPlaque_filter$HGNC %in% geneSub])
map4 <- match(geneSub, resPlaque_cells_filter$HGNC[resPlaque_cells_filter$HGNC %in% geneSub])

baseMeans <- c(resPlaqueInd_filter$baseMean[resPlaqueInd_filter$HGNC %in% geneSub][map1],
               resPlaqueInd_cells_filter$baseMean[resPlaqueInd_cells_filter$HGNC %in% geneSub][map2],
               resPlaque_filter$baseMean[resPlaque_filter$HGNC %in% geneSub][map3],
               resPlaque_cells_filter$baseMean[resPlaque_cells_filter$HGNC %in% geneSub][map4])

l2fc <- c(resPlaqueInd_filter$log2FoldChange[resPlaqueInd_filter$HGNC %in% geneSub][map1],
               resPlaqueInd_cells_filter$log2FoldChange[resPlaqueInd_cells_filter$HGNC %in% geneSub][map2],
               resPlaque_filter$log2FoldChange[resPlaque_filter$HGNC %in% geneSub][map3],
               resPlaque_cells_filter$log2FoldChange[resPlaque_cells_filter$HGNC %in% geneSub][map4]) %>% round(digits = 4)

l2fcse <- c(resPlaqueInd_filter$lfcSE[resPlaqueInd_filter$HGNC %in% geneSub][map1],
               resPlaqueInd_cells_filter$lfcSE[resPlaqueInd_cells_filter$HGNC %in% geneSub][map2],
               resPlaque_filter$lfcSE[resPlaque_filter$HGNC %in% geneSub][map3],
               resPlaque_cells_filter$lfcSE[resPlaque_cells_filter$HGNC %in% geneSub][map4]) %>% round(digits = 4)

wStat <- c(resPlaqueInd_filter$stat[resPlaqueInd_filter$HGNC %in% geneSub][map1],
            resPlaqueInd_cells_filter$stat[resPlaqueInd_cells_filter$HGNC %in% geneSub][map2],
            resPlaque_filter$stat[resPlaque_filter$HGNC %in% geneSub][map3],
            resPlaque_cells_filter$stat[resPlaque_cells_filter$HGNC %in% geneSub][map4]) %>% round(digits = 4)

padjs <- c(resPlaqueInd_filter$padj[resPlaqueInd_filter$HGNC %in% geneSub][map1],
               resPlaqueInd_cells_filter$padj[resPlaqueInd_cells_filter$HGNC %in% geneSub][map2],
               resPlaque_filter$padj[resPlaque_filter$HGNC %in% geneSub][map3],
               resPlaque_cells_filter$padj[resPlaque_cells_filter$HGNC %in% geneSub][map4]) %>% round(digits = 4)

plaqueDf <- data.frame(Gene = genesDouble,
                       Association = asso,
                       Adjusted = adjusted,
                       Plaques = plaque,
                       BaseMean = baseMeans,
                       log2FoldChange = l2fc,
                       log2FoldChangeSE = l2fcse,
                       stat = wStat,
                       pAdjust = padjs) %>% arrange(plaqueDf, Gene, Adjusted, Plaques)

##########################################
####################3##################
map1 <- match(geneSub, resPlaqueInd$HGNC[resPlaqueInd$HGNC %in% geneSub])
map2 <- match(geneSub, resPlaqueInd_cells$HGNC[resPlaqueInd_cells$HGNC %in% geneSub])
map3 <- match(geneSub, resPlaque$HGNC[resPlaque$HGNC %in% geneSub])
map4 <- match(geneSub, resPlaque_cells$HGNC[resPlaque_cells$HGNC %in% geneSub])

baseMeans <- c(resPlaqueInd$baseMean[resPlaqueInd$HGNC %in% geneSub][map1],
               resPlaqueInd_cells$baseMean[resPlaqueInd_cells$HGNC %in% geneSub][map2],
               resPlaque$baseMean[resPlaque$HGNC %in% geneSub][map3],
               resPlaque_cells$baseMean[resPlaque_cells$HGNC %in% geneSub][map4])

l2fc <- c(resPlaqueInd$log2FoldChange[resPlaqueInd$HGNC %in% geneSub][map1],
          resPlaqueInd_cells$log2FoldChange[resPlaqueInd_cells$HGNC %in% geneSub][map2],
          resPlaque$log2FoldChange[resPlaque$HGNC %in% geneSub][map3],
          resPlaque_cells$log2FoldChange[resPlaque_cells$HGNC %in% geneSub][map4]) %>% round(digits = 4)

l2fcse <- c(resPlaqueInd$lfcSE[resPlaqueInd$HGNC %in% geneSub][map1],
            resPlaqueInd_cells$lfcSE[resPlaqueInd_cells$HGNC %in% geneSub][map2],
            resPlaque$lfcSE[resPlaque$HGNC %in% geneSub][map3],
            resPlaque_cells$lfcSE[resPlaque_cells$HGNC %in% geneSub][map4]) %>% round(digits = 4)

wStat <- c(resPlaqueInd$stat[resPlaqueInd$HGNC %in% geneSub][map1],
           resPlaqueInd_cells$stat[resPlaqueInd_cells$HGNC %in% geneSub][map2],
           resPlaque$stat[resPlaque$HGNC %in% geneSub][map3],
           resPlaque_cells$stat[resPlaque_cells$HGNC %in% geneSub][map4]) %>% round(digits = 4)

padjs <- c(resPlaqueInd$padj[resPlaqueInd$HGNC %in% geneSub][map1],
           resPlaqueInd_cells$padj[resPlaqueInd_cells$HGNC %in% geneSub][map2],
           resPlaque$padj[resPlaque$HGNC %in% geneSub][map3],
           resPlaque_cells$padj[resPlaque_cells$HGNC %in% geneSub][map4]) %>% round(digits = 4)

plaqueDf2 <- data.frame(Gene = genesDouble,
                       #Association = asso,
                       Adjusted = adjusted,
                       Plaques = plaque,
                       BaseMean = baseMeans,
                       log2FoldChange = l2fc,
                       log2FoldChangeSE = l2fcse,
                       stat = wStat,
                       pAdjust = padjs) %>% arrange(., Gene, Adjusted, Plaques)

##########################################
##########################################

duplicatedGenes <- resPlaque_filter$HGNC[which(duplicated(resPlaque_filter$HGNC))]

resPlaque_filter$HGNC[which(resPlaque_filter$HGNC == "SCO2")] <- c("SCO2-1", "SCO2-2")

resPlaque_filter$HGNC[which(resPlaque_filter$HGNC == "RABGEF1")] <- c("RABGEF1-1", "RABGEF1-2")

resPlaque_filter$HGNC[which(resPlaque_filter$HGNC == "MATR3")] <- c("MATR3-1", "MATR3-2")

resPlaque_filter$HGNC[which(resPlaque_filter$HGNC == "MALAT1")] <- c("MALAT1-1", "MALAT2-2")

resPlaque_filter$HGNC[which(resPlaque_filter$HGNC == "COG8")] <- c("COG8-1", "COG8-2")

resPlaque_cells_filter$HGNC[which(resPlaque_cells_filter$HGNC == "SCO2")] <- c("SCO2-1", "SCO2-2")

resPlaque_cells_filter$HGNC[which(resPlaque_cells_filter$HGNC == "RABGEF1")] <- c("RABGEF1-1", "RABGEF1-2")

resPlaque_cells_filter$HGNC[which(resPlaque_cells_filter$HGNC == "MATR3")] <- c("MATR3-1", "MATR3-2")

resPlaque_cells_filter$HGNC[which(resPlaque_cells_filter$HGNC == "MALAT1")] <- c("MALAT1-1", "MALAT2-2")

resPlaque_cells_filter$HGNC[which(resPlaque_cells_filter$HGNC == "COG8")] <- c("COG8-1", "COG8-2")

up_up_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$log2FoldChange > 0 & resPlaque_filter$padj <= 0.1], 
                         resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$log2FoldChange > 0 & resPlaque_cells_filter$padj <= 0.1])

up_down_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$log2FoldChange > 0 & resPlaque_filter$padj <= 0.1], 
                           resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$log2FoldChange < 0 & resPlaque_cells_filter$padj <= 0.1])

down_down_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$log2FoldChange < 0 & resPlaque_filter$padj <= 0.1], 
                           resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$log2FoldChange < 0 & resPlaque_cells_filter$padj <= 0.1])

down_up_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$log2FoldChange < 0 & resPlaque_filter$padj <= 0.1], 
                           resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$log2FoldChange > 0 & resPlaque_cells_filter$padj <= 0.1])

up_no_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$log2FoldChange > 0 & resPlaque_filter$padj <= 0.1], 
                           resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$padj > 0.1])

down_no_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$log2FoldChange < 0 & resPlaque_filter$padj <= 0.1], 
                         resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$padj > 0.1])

no_up_genes <- intersect(resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$log2FoldChange > 0 & resPlaque_cells_filter$padj <= 0.1], 
                         resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$padj > 0.1])

no_down_genes <- intersect(resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$log2FoldChange < 0 & resPlaque_cells_filter$padj <= 0.1], 
                         resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$padj > 0.1])

no_no_genes <- intersect(resPlaque_filter$HGNC[!(is.na(resPlaque_filter$padj)) & resPlaque_filter$padj > 0.1],
                         resPlaque_cells_filter$HGNC[!(is.na(resPlaque_cells_filter$padj)) & resPlaque_cells_filter$padj > 0.1])

up_up_df <- data.frame(Gene = up_up_genes, 
                       BaseMean = vector(mode = "double", length = length(up_up_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(up_up_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(up_up_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(up_up_genes)),
                       AdjustedStat = vector(mode = "double", length = length(up_up_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(up_up_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(up_up_genes)))
for(i in 1:length(up_up_genes)){
  ind1 <- which(resPlaque_filter$HGNC == up_up_genes[i])
  ind2 <- which(resPlaque_cells_filter$HGNC == up_up_genes[i])
  
  up_up_df$BaseMean[i] <- resPlaque_filter$baseMean[ind1] %>% round(digits = 4)
  up_up_df$UnadjustedLFC[i] <- resPlaque_filter$log2FoldChange[ind1] %>% round(digits = 4)
  up_up_df$AdjustedLFC[i] <- resPlaque_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  up_up_df$UnadjustedStat[i] <- resPlaque_filter$stat[ind1] %>% round(digits = 4)
  up_up_df$AdjustedStat[i] <- resPlaque_cells_filter$stat[ind2] %>% round(digits = 4)
  up_up_df$UnadjustedPadj[i] <- resPlaque_filter$padj[ind1] %>% round(digits = 4)
  up_up_df$AdjustedPadj[i] <- resPlaque_cells_filter$padj[ind2] %>% round(digits = 4)
}

down_down_df <- data.frame(Gene = down_down_genes, 
                       BaseMean = vector(mode = "double", length = length(down_down_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(down_down_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(down_down_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(down_down_genes)),
                       AdjustedStat = vector(mode = "double", length = length(down_down_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(down_down_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(down_down_genes)))
for(i in 1:length(down_down_genes)){
  ind1 <- which(resPlaque_filter$HGNC == down_down_genes[i])
  ind2 <- which(resPlaque_cells_filter$HGNC == down_down_genes[i])
  
  down_down_df$BaseMean[i] <- resPlaque_filter$baseMean[ind1] %>% round(digits = 4)
  down_down_df$UnadjustedLFC[i] <- resPlaque_filter$log2FoldChange[ind1] %>% round(digits = 4)
  down_down_df$AdjustedLFC[i] <- resPlaque_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  down_down_df$UnadjustedStat[i] <- resPlaque_filter$stat[ind1] %>% round(digits = 4)
  down_down_df$AdjustedStat[i] <- resPlaque_cells_filter$stat[ind2] %>% round(digits = 4)
  down_down_df$UnadjustedPadj[i] <- resPlaque_filter$padj[ind1] %>% round(digits = 4)
  down_down_df$AdjustedPadj[i] <- resPlaque_cells_filter$padj[ind2] %>% round(digits = 4)
}

no_up_df <- data.frame(Gene = no_up_genes, 
                       BaseMean = vector(mode = "double", length = length(no_up_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(no_up_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(no_up_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(no_up_genes)),
                       AdjustedStat = vector(mode = "double", length = length(no_up_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(no_up_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(no_up_genes)))
for(i in 1:length(no_up_genes)){
  ind1 <- which(resPlaque_filter$HGNC == no_up_genes[i])
  ind2 <- which(resPlaque_cells_filter$HGNC == no_up_genes[i])
  
  no_up_df$BaseMean[i] <- resPlaque_filter$baseMean[ind1] %>% round(digits = 4)
  no_up_df$UnadjustedLFC[i] <- resPlaque_filter$log2FoldChange[ind1] %>% round(digits = 4)
  no_up_df$AdjustedLFC[i] <- resPlaque_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  no_up_df$UnadjustedStat[i] <- resPlaque_filter$stat[ind1] %>% round(digits = 4)
  no_up_df$AdjustedStat[i] <- resPlaque_cells_filter$stat[ind2] %>% round(digits = 4)
  no_up_df$UnadjustedPadj[i] <- resPlaque_filter$padj[ind1] %>% round(digits = 4)
  no_up_df$AdjustedPadj[i] <- resPlaque_cells_filter$padj[ind2] %>% round(digits = 4)
}

no_down_df <- data.frame(Gene = no_down_genes, 
                       BaseMean = vector(mode = "double", length = length(no_down_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(no_down_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(no_down_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(no_down_genes)),
                       AdjustedStat = vector(mode = "double", length = length(no_down_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(no_down_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(no_down_genes)))
for(i in 1:length(no_down_genes)){
  ind1 <- which(resPlaque_filter$HGNC == no_down_genes[i])
  ind2 <- which(resPlaque_cells_filter$HGNC == no_down_genes[i])
  
  no_down_df$BaseMean[i] <- resPlaque_filter$baseMean[ind1] %>% round(digits = 4)
  no_down_df$UnadjustedLFC[i] <- resPlaque_filter$log2FoldChange[ind1] %>% round(digits = 4)
  no_down_df$AdjustedLFC[i] <- resPlaque_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  no_down_df$UnadjustedStat[i] <- resPlaque_filter$stat[ind1] %>% round(digits = 4)
  no_down_df$AdjustedStat[i] <- resPlaque_cells_filter$stat[ind2] %>% round(digits = 4)
  no_down_df$UnadjustedPadj[i] <- resPlaque_filter$padj[ind1] %>% round(digits = 4)
  no_down_df$AdjustedPadj[i] <- resPlaque_cells_filter$padj[ind2] %>% round(digits = 4)
}

up_no_df <- data.frame(Gene = up_no_genes, 
                       BaseMean = vector(mode = "double", length = length(up_no_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(up_no_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(up_no_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(up_no_genes)),
                       AdjustedStat = vector(mode = "double", length = length(up_no_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(up_no_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(up_no_genes)))
for(i in 1:length(up_no_genes)){
  ind1 <- which(resPlaque_filter$HGNC == up_no_genes[i])
  ind2 <- which(resPlaque_cells_filter$HGNC == up_no_genes[i])
  
  up_no_df$BaseMean[i] <- resPlaque_filter$baseMean[ind1] %>% round(digits = 4)
  up_no_df$UnadjustedLFC[i] <- resPlaque_filter$log2FoldChange[ind1] %>% round(digits = 4)
  up_no_df$AdjustedLFC[i] <- resPlaque_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  up_no_df$UnadjustedStat[i] <- resPlaque_filter$stat[ind1] %>% round(digits = 4)
  up_no_df$AdjustedStat[i] <- resPlaque_cells_filter$stat[ind2] %>% round(digits = 4)
  up_no_df$UnadjustedPadj[i] <- resPlaque_filter$padj[ind1] %>% round(digits = 4)
  up_no_df$AdjustedPadj[i] <- resPlaque_cells_filter$padj[ind2] %>% round(digits = 4)
}

down_no_df <- data.frame(Gene = down_no_genes, 
                       BaseMean = vector(mode = "double", length = length(down_no_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(down_no_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(down_no_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(down_no_genes)),
                       AdjustedStat = vector(mode = "double", length = length(down_no_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(down_no_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(down_no_genes)))
for(i in 1:length(down_no_genes)){
  ind1 <- which(resPlaque_filter$HGNC == down_no_genes[i])
  ind2 <- which(resPlaque_cells_filter$HGNC == down_no_genes[i])
  
  down_no_df$BaseMean[i] <- resPlaque_filter$baseMean[ind1] %>% round(digits = 4)
  down_no_df$UnadjustedLFC[i] <- resPlaque_filter$log2FoldChange[ind1] %>% round(digits = 4)
  down_no_df$AdjustedLFC[i] <- resPlaque_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  down_no_df$UnadjustedStat[i] <- resPlaque_filter$stat[ind1] %>% round(digits = 4)
  down_no_df$AdjustedStat[i] <- resPlaque_cells_filter$stat[ind2] %>% round(digits = 4)
  down_no_df$UnadjustedPadj[i] <- resPlaque_filter$padj[ind1] %>% round(digits = 4)
  down_no_df$AdjustedPadj[i] <- resPlaque_cells_filter$padj[ind2] %>% round(digits = 4)
}

############################################
############################################

##########################################
##########################################

duplicatedGenes <- resPlaqueInd_filter$HGNC[which(duplicated(resPlaqueInd_filter$HGNC))]

resPlaqueInd_filter$HGNC[which(resPlaqueInd_filter$HGNC == "SCO2")] <- c("SCO2-1", "SCO2-2")

resPlaqueInd_filter$HGNC[which(resPlaqueInd_filter$HGNC == "RABGEF1")] <- c("RABGEF1-1", "RABGEF1-2")

resPlaqueInd_filter$HGNC[which(resPlaqueInd_filter$HGNC == "MATR3")] <- c("MATR3-1", "MATR3-2")

resPlaqueInd_filter$HGNC[which(resPlaqueInd_filter$HGNC == "MALAT1")] <- c("MALAT1-1", "MALAT2-2")

resPlaqueInd_filter$HGNC[which(resPlaqueInd_filter$HGNC == "COG8")] <- c("COG8-1", "COG8-2")

resPlaqueInd_cells_filter$HGNC[which(resPlaqueInd_cells_filter$HGNC == "SCO2")] <- c("SCO2-1", "SCO2-2")

resPlaqueInd_cells_filter$HGNC[which(resPlaqueInd_cells_filter$HGNC == "RABGEF1")] <- c("RABGEF1-1", "RABGEF1-2")

resPlaqueInd_cells_filter$HGNC[which(resPlaqueInd_cells_filter$HGNC == "MATR3")] <- c("MATR3-1", "MATR3-2")

resPlaqueInd_cells_filter$HGNC[which(resPlaqueInd_cells_filter$HGNC == "MALAT1")] <- c("MALAT1-1", "MALAT2-2")

resPlaqueInd_cells_filter$HGNC[which(resPlaqueInd_cells_filter$HGNC == "COG8")] <- c("COG8-1", "COG8-2")

up_up_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$log2FoldChange > 0 & resPlaqueInd_filter$padj <= 0.1], 
                         resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$log2FoldChange > 0 & resPlaqueInd_cells_filter$padj <= 0.1])

up_down_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$log2FoldChange > 0 & resPlaqueInd_filter$padj <= 0.1], 
                           resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$log2FoldChange < 0 & resPlaqueInd_cells_filter$padj <= 0.1])

down_down_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$log2FoldChange < 0 & resPlaqueInd_filter$padj <= 0.1], 
                             resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$log2FoldChange < 0 & resPlaqueInd_cells_filter$padj <= 0.1])

down_up_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$log2FoldChange < 0 & resPlaqueInd_filter$padj <= 0.1], 
                           resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$log2FoldChange > 0 & resPlaqueInd_cells_filter$padj <= 0.1])

up_no_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$log2FoldChange > 0 & resPlaqueInd_filter$padj <= 0.1], 
                         resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$padj > 0.1])

down_no_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$log2FoldChange < 0 & resPlaqueInd_filter$padj <= 0.1], 
                           resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$padj > 0.1])

no_up_genes <- intersect(resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$log2FoldChange > 0 & resPlaqueInd_cells_filter$padj <= 0.1], 
                         resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$padj > 0.1])

no_down_genes <- intersect(resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$log2FoldChange < 0 & resPlaqueInd_cells_filter$padj <= 0.1], 
                           resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$padj > 0.1])

no_no_genes <- intersect(resPlaqueInd_filter$HGNC[!(is.na(resPlaqueInd_filter$padj)) & resPlaqueInd_filter$padj > 0.1],
                         resPlaqueInd_cells_filter$HGNC[!(is.na(resPlaqueInd_cells_filter$padj)) & resPlaqueInd_cells_filter$padj > 0.1])

up_up_df <- data.frame(Gene = up_up_genes, 
                       BaseMean = vector(mode = "double", length = length(up_up_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(up_up_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(up_up_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(up_up_genes)),
                       AdjustedStat = vector(mode = "double", length = length(up_up_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(up_up_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(up_up_genes)))
for(i in 1:length(up_up_genes)){
  ind1 <- which(resPlaqueInd_filter$HGNC == up_up_genes[i])
  ind2 <- which(resPlaqueInd_cells_filter$HGNC == up_up_genes[i])
  
  up_up_df$BaseMean[i] <- resPlaqueInd_filter$baseMean[ind1] %>% round(digits = 4)
  up_up_df$UnadjustedLFC[i] <- resPlaqueInd_filter$log2FoldChange[ind1] %>% round(digits = 4)
  up_up_df$AdjustedLFC[i] <- resPlaqueInd_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  up_up_df$UnadjustedStat[i] <- resPlaqueInd_filter$stat[ind1] %>% round(digits = 4)
  up_up_df$AdjustedStat[i] <- resPlaqueInd_cells_filter$stat[ind2] %>% round(digits = 4)
  up_up_df$UnadjustedPadj[i] <- resPlaqueInd_filter$padj[ind1] %>% round(digits = 4)
  up_up_df$AdjustedPadj[i] <- resPlaqueInd_cells_filter$padj[ind2] %>% round(digits = 4)
}

down_down_df <- data.frame(Gene = down_down_genes, 
                           BaseMean = vector(mode = "double", length = length(down_down_genes)),
                           UnadjustedLFC = vector(mode = "double", length = length(down_down_genes)),
                           AdjustedLFC = vector(mode = "double", length = length(down_down_genes)),
                           UnadjustedStat = vector(mode = "double", length = length(down_down_genes)),
                           AdjustedStat = vector(mode = "double", length = length(down_down_genes)),
                           UnadjustedPadj = vector(mode = "double", length = length(down_down_genes)),
                           AdjustedPadj = vector(mode = "double", length = length(down_down_genes)))
for(i in 1:length(down_down_genes)){
  ind1 <- which(resPlaqueInd_filter$HGNC == down_down_genes[i])
  ind2 <- which(resPlaqueInd_cells_filter$HGNC == down_down_genes[i])
  
  down_down_df$BaseMean[i] <- resPlaqueInd_filter$baseMean[ind1] %>% round(digits = 4)
  down_down_df$UnadjustedLFC[i] <- resPlaqueInd_filter$log2FoldChange[ind1] %>% round(digits = 4)
  down_down_df$AdjustedLFC[i] <- resPlaqueInd_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  down_down_df$UnadjustedStat[i] <- resPlaqueInd_filter$stat[ind1] %>% round(digits = 4)
  down_down_df$AdjustedStat[i] <- resPlaqueInd_cells_filter$stat[ind2] %>% round(digits = 4)
  down_down_df$UnadjustedPadj[i] <- resPlaqueInd_filter$padj[ind1] %>% round(digits = 4)
  down_down_df$AdjustedPadj[i] <- resPlaqueInd_cells_filter$padj[ind2] %>% round(digits = 4)
}

no_up_df <- data.frame(Gene = no_up_genes, 
                       BaseMean = vector(mode = "double", length = length(no_up_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(no_up_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(no_up_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(no_up_genes)),
                       AdjustedStat = vector(mode = "double", length = length(no_up_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(no_up_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(no_up_genes)))
for(i in 1:length(no_up_genes)){
  ind1 <- which(resPlaqueInd_filter$HGNC == no_up_genes[i])
  ind2 <- which(resPlaqueInd_cells_filter$HGNC == no_up_genes[i])
  
  no_up_df$BaseMean[i] <- resPlaqueInd_filter$baseMean[ind1] %>% round(digits = 4)
  no_up_df$UnadjustedLFC[i] <- resPlaqueInd_filter$log2FoldChange[ind1] %>% round(digits = 4)
  no_up_df$AdjustedLFC[i] <- resPlaqueInd_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  no_up_df$UnadjustedStat[i] <- resPlaqueInd_filter$stat[ind1] %>% round(digits = 4)
  no_up_df$AdjustedStat[i] <- resPlaqueInd_cells_filter$stat[ind2] %>% round(digits = 4)
  no_up_df$UnadjustedPadj[i] <- resPlaqueInd_filter$padj[ind1] %>% round(digits = 4)
  no_up_df$AdjustedPadj[i] <- resPlaqueInd_cells_filter$padj[ind2] %>% round(digits = 4)
}

no_down_df <- data.frame(Gene = no_down_genes, 
                         BaseMean = vector(mode = "double", length = length(no_down_genes)),
                         UnadjustedLFC = vector(mode = "double", length = length(no_down_genes)),
                         AdjustedLFC = vector(mode = "double", length = length(no_down_genes)),
                         UnadjustedStat = vector(mode = "double", length = length(no_down_genes)),
                         AdjustedStat = vector(mode = "double", length = length(no_down_genes)),
                         UnadjustedPadj = vector(mode = "double", length = length(no_down_genes)),
                         AdjustedPadj = vector(mode = "double", length = length(no_down_genes)))
for(i in 1:length(no_down_genes)){
  ind1 <- which(resPlaqueInd_filter$HGNC == no_down_genes[i])
  ind2 <- which(resPlaqueInd_cells_filter$HGNC == no_down_genes[i])
  
  no_down_df$BaseMean[i] <- resPlaqueInd_filter$baseMean[ind1] %>% round(digits = 4)
  no_down_df$UnadjustedLFC[i] <- resPlaqueInd_filter$log2FoldChange[ind1] %>% round(digits = 4)
  no_down_df$AdjustedLFC[i] <- resPlaqueInd_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  no_down_df$UnadjustedStat[i] <- resPlaqueInd_filter$stat[ind1] %>% round(digits = 4)
  no_down_df$AdjustedStat[i] <- resPlaqueInd_cells_filter$stat[ind2] %>% round(digits = 4)
  no_down_df$UnadjustedPadj[i] <- resPlaqueInd_filter$padj[ind1] %>% round(digits = 4)
  no_down_df$AdjustedPadj[i] <- resPlaqueInd_cells_filter$padj[ind2] %>% round(digits = 4)
}

up_no_df <- data.frame(Gene = up_no_genes, 
                       BaseMean = vector(mode = "double", length = length(up_no_genes)),
                       UnadjustedLFC = vector(mode = "double", length = length(up_no_genes)),
                       AdjustedLFC = vector(mode = "double", length = length(up_no_genes)),
                       UnadjustedStat = vector(mode = "double", length = length(up_no_genes)),
                       AdjustedStat = vector(mode = "double", length = length(up_no_genes)),
                       UnadjustedPadj = vector(mode = "double", length = length(up_no_genes)),
                       AdjustedPadj = vector(mode = "double", length = length(up_no_genes)))
for(i in 1:length(up_no_genes)){
  ind1 <- which(resPlaqueInd_filter$HGNC == up_no_genes[i])
  ind2 <- which(resPlaqueInd_cells_filter$HGNC == up_no_genes[i])
  
  up_no_df$BaseMean[i] <- resPlaqueInd_filter$baseMean[ind1] %>% round(digits = 4)
  up_no_df$UnadjustedLFC[i] <- resPlaqueInd_filter$log2FoldChange[ind1] %>% round(digits = 4)
  up_no_df$AdjustedLFC[i] <- resPlaqueInd_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  up_no_df$UnadjustedStat[i] <- resPlaqueInd_filter$stat[ind1] %>% round(digits = 4)
  up_no_df$AdjustedStat[i] <- resPlaqueInd_cells_filter$stat[ind2] %>% round(digits = 4)
  up_no_df$UnadjustedPadj[i] <- resPlaqueInd_filter$padj[ind1] %>% round(digits = 4)
  up_no_df$AdjustedPadj[i] <- resPlaqueInd_cells_filter$padj[ind2] %>% round(digits = 4)
}

down_no_df <- data.frame(Gene = down_no_genes, 
                         BaseMean = vector(mode = "double", length = length(down_no_genes)),
                         UnadjustedLFC = vector(mode = "double", length = length(down_no_genes)),
                         AdjustedLFC = vector(mode = "double", length = length(down_no_genes)),
                         UnadjustedStat = vector(mode = "double", length = length(down_no_genes)),
                         AdjustedStat = vector(mode = "double", length = length(down_no_genes)),
                         UnadjustedPadj = vector(mode = "double", length = length(down_no_genes)),
                         AdjustedPadj = vector(mode = "double", length = length(down_no_genes)))
for(i in 1:length(down_no_genes)){
  ind1 <- which(resPlaqueInd_filter$HGNC == down_no_genes[i])
  ind2 <- which(resPlaqueInd_cells_filter$HGNC == down_no_genes[i])
  
  down_no_df$BaseMean[i] <- resPlaqueInd_filter$baseMean[ind1] %>% round(digits = 4)
  down_no_df$UnadjustedLFC[i] <- resPlaqueInd_filter$log2FoldChange[ind1] %>% round(digits = 4)
  down_no_df$AdjustedLFC[i] <- resPlaqueInd_cells_filter$log2FoldChange[ind2] %>% round(digits = 4)
  down_no_df$UnadjustedStat[i] <- resPlaqueInd_filter$stat[ind1] %>% round(digits = 4)
  down_no_df$AdjustedStat[i] <- resPlaqueInd_cells_filter$stat[ind2] %>% round(digits = 4)
  down_no_df$UnadjustedPadj[i] <- resPlaqueInd_filter$padj[ind1] %>% round(digits = 4)
  down_no_df$AdjustedPadj[i] <- resPlaqueInd_cells_filter$padj[ind2] %>% round(digits = 4)
}


####################################
####################################

plaqueRes <- resPlaque_filter %>% data.frame()
colnames(plaqueRes) <- c("HGNC", "baseMean", "l2fc_U", "l2fcSE_U", "stat_U", "pval_U", "padj_U")
map <- match(plaqueRes$HGNC, resPlaque_cells_filter$HGNC[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC])
identical(resPlaque_cells_filter$HGNC[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC][map], plaqueRes$HGNC)
plaqueRes <- mutate(plaqueRes, 
                    l2fc_A = resPlaque_cells_filter$log2FoldChange[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC][map],
                    l2fcSE_A = resPlaque_cells_filter$lfcSE[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC][map],
                    stat_A = resPlaque_cells_filter$stat[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC][map],
                    pval_A = resPlaque_cells_filter$pvalue[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC][map],
                    padj_A = resPlaque_cells_filter$padj[resPlaque_cells_filter$HGNC %in% plaqueRes$HGNC][map],
                    l2fcDiff = l2fc_A - l2fc_U) %>% arrange(baseMean)

pA_minus_U <- ggplot(plaqueRes, aes(x = log2(baseMean), y = l2fcDiff, text = HGNC, color = abs(l2fcDiff) >= 0.5, alpha = 0.5)) + geom_point() +
                labs(title = "Plaque 4 vs 0") + scale_color_manual(values = c("black", "red")) + 
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", text = element_text(family = "Fira Code"))

####################################
####################################

plaqueResInd <- resPlaqueInd_filter %>% data.frame()
colnames(plaqueResInd) <- c("HGNC", "baseMean", "l2fc_U", "l2fcSE_U", "stat_U", "pval_U", "padj_U")
map <- match(plaqueResInd$HGNC, resPlaqueInd_cells_filter$HGNC[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC])
identical(resPlaqueInd_cells_filter$HGNC[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC][map], plaqueResInd$HGNC)
plaqueResInd <- mutate(plaqueResInd, 
                    l2fc_A = resPlaqueInd_cells_filter$log2FoldChange[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC][map],
                    l2fcSE_A = resPlaqueInd_cells_filter$lfcSE[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC][map],
                    stat_A = resPlaqueInd_cells_filter$stat[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC][map],
                    pval_A = resPlaqueInd_cells_filter$pvalue[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC][map],
                    padj_A = resPlaqueInd_cells_filter$padj[resPlaqueInd_cells_filter$HGNC %in% plaqueResInd$HGNC][map],
                    l2fcDiff = l2fc_A - l2fc_U) %>% arrange(baseMean)

pA_minus_UInd <- ggplot(plaqueResInd, aes(x = log2(baseMean), y = l2fcDiff, text = HGNC, color = abs(l2fcDiff) >= 0.5, alpha = 0.5)) + geom_point() +
  labs(title = "All vs 0", y = "l2FC_A - l2FC_U") + scale_color_manual(values = c("black", "red")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", text = element_text(family = "Fira Code"))

ggplotly(pA_minus_U, tooltip = c("x", "y", "text"))

ggplotly(pA_minus_UInd, tooltip = c("x", "y", "text"))

#################################
#################################


q_ij <- t( t( assays(deseqPlaque)[["mu"]] ) / sizeFactors(deseqPlaque) )
mod <- attr(deseqPlaque, "modelMatrix")
coefs <- coef(deseqPlaque)
q_tst <- coefs[,-11] %*% t(mod)[-11,]

df <- data.frame(sample = colnames(deseqPlaque), log2q = log2(q_ij)[6474,], qtest = q_tst[6474,],
                 GeneExp = counts(deseqPlaque_cells, normalized = TRUE)[6474,], expQ = exp(q_tst_cells[6474,]),
                 plaque = factor(deseqPlaque$Plaque == 4, labels = c(0, 4)),
                 end = deseqPlaque_cells$END, mac = deseqPlaque_cells$MAC) %>% arrange(mac)

ggplot(df, aes(x = sample, y = GeneExp - expQ, color = mac)) + geom_point() + labs(title = "APOC1", x = "Sample", y = "Adjusted Expression") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + scale_color_viridis_c() 




q_ij_cells <- t( t( assays(deseqPlaque_cells)[["mu"]] ) / sizeFactors(deseqPlaque_cells) )
mod_cells <- attr(deseqPlaque_cells, "modelMatrix")
coefs_cells <- coef(deseqPlaque_cells)
q_tst_cells <- coefs_cells[,-5] %*% t(mod_cells)[-5,]

df_cells <- data.frame(sample = colnames(deseqPlaque_cells), log2q = log2(q_ij_cells)[6474,], qtest = q_tst_cells[6474,], plaque = factor(deseqPlaque_cells$Plaque == 4))

ggplot(df_cells, aes(x = sample, y = log2q - qtest, color = plaque)) + geom_point() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + scale_color_manual(values = c("black", "red")) + 















############################################
############################################

plaqueRes <- results(gtexPlaque, contrast = c("Plaque", 4, 0))
plaqueRes$HGNC <- as.character(df)
plaqueRes <- plaqueRes[,c(7, 1:6)]
plaqueResSort <- plaqueRes[order(plaqueRes$padj, decreasing = FALSE), ]

pcellsRes <- results(gtexPCells, contrast = c("Plaque", 4, 0))
pcellsRes$HGNC <- as.character(df)
pcellsRes <- pcellsRes[,c(7, 1:6)]
pcellsSort <- pcellsRes[order(pcellsRes$padj, decreasing = FALSE), ]


proportionFactors2 <- proportionFactors1
proportionFactors2$ADI <- (celltypeProportions$ADI > cellPropSummaries["3rd Qu.", "ADI"])*1
gtex$ADI <- factor(proportionFactors2$ADI)
gtexPCells2 <- DESeqDataSet(gtex, design = ~ Plaque + CMC + END + LYM + SMC + MAC + ADI) %>% DESeq()
pcellsRes2 <- results(gtexPCells2, contrast = c("Plaque", 4, 0))
pcellsRes2$HGNC <- as.character(df)
pcellsRes2 <- pcellsRes2[,c(7, 1:6)]
pcellsSort2 <- pcellsRes2[order(pcellsRes2$padj, decreasing = FALSE), ]

proportionFactors3 <- proportionFactors1
proportionFactors3$ADI <- (celltypeProportions$ADI > cellPropSummaries["3rd Qu.", "ADI"])*1
gtex$ADI <- factor(proportionFactors2$ADI)
gtexPCells2 <- DESeqDataSet(gtex, design = ~ Plaque + CMC + END + LYM + SMC + MAC + ADI) %>% DESeq()
pcellsRes2 <- results(gtexPCells2, contrast = c("Plaque", 4, 0))
pcellsRes2$HGNC <- as.character(df)
pcellsRes2 <- pcellsRes2[,c(7, 1:6)]
pcellsSort2 <- pcellsRes2[order(pcellsRes2$padj, decreasing = FALSE), ]

gtexPCells4 <- DESeqDataSet(gtex, design = ~ Plaque + CMC + ADI) %>% DESeq()
pcellsRes4 <- results(gtexPCells4, contrast = c("Plaque", 4, 0))
pcellsRes4$HGNC <- as.character(df)
pcellsRes4 <- pcellsRes4[,c(7, 1:6)]
pcellsSort4 <- pcellsRes4[order(pcellsRes4$padj, decreasing = FALSE), ]

cmc_adi <- celltypeProportions$CMC + celltypeProportions$ADI
end_smc <- celltypeProportions$SMC + celltypeProportions$END
lym_mac <- celltypeProportions$LYM

inf <- as.integer(lym_mac > summary(lym_mac)["1st Qu."])
inf[which(lym_mac > summary(lym_mac)["3rd Qu."])] <- 2
ext <- as.integer(cmc_adi > summary(cmc_adi)["Median"])
ext[which(cmc_adi > summary(cmc_adi)["3rd Qu."])] <- 2
nml <- as.integer(end_smc > summary(end_smc)["Median"])
nml[which(end_smc > summary(end_smc)["3rd Qu."])] <- 2
gtex$EXT <- factor(ext)
gtex$INF <- factor(inf)
gtex$NML <- factor(nml)
gtexPCells4 <- DESeqDataSet(gtex, design = ~ Plaque + INF + EXT) %>% DESeq()
pcellsRes4 <- results(gtexPCells4, contrast = c("Plaque", 4, 0))
pcellsRes4$HGNC <- as.character(df)
pcellsRes4 <- pcellsRes4[,c(7, 1:6)]
pcellsSort4 <- pcellsRes4[order(pcellsRes4$padj, decreasing = FALSE), ]

gtexPCells5 <- DESeqDataSet(gtex, design = ~ Plaque + NML + EXT + INF) %>% DESeq()
pcellsRes5 <- results(gtexPCells5, contrast = c("Plaque", 4, 0))
pcellsRes5$HGNC <- as.character(df)
pcellsRes5 <- pcellsRes5[,c(7, 1:6)]
pcellsSort5 <- pcellsRes5[order(pcellsRes5$padj, decreasing = FALSE), ]

gtexPCells6 <- DESeqDataSet(gtex, design = ~ Plaque + EXT + SMC) %>% DESeq()
pcellsRes6 <- results(gtexPCells6, contrast = c("Plaque", 4, 0))
pcellsRes6$HGNC <- as.character(df)
pcellsRes6 <- pcellsRes6[,c(7, 1:6)]
pcellsSort6 <- pcellsRes6[order(pcellsRes6$padj, decreasing = FALSE), ]

gtexPCells7 <- DESeqDataSet(gtex, design = ~ Plaque + CMC + END + LYM + SMC + ADI) %>% DESeq()
pcellsRes7 <- results(gtexPCells7, contrast = c("Plaque", 4, 0))
pcellsRes7$HGNC <- as.character(df)
pcellsRes7 <- pcellsRes7[,c(7, 1:6)]
pcellsSort7 <- pcellsRes7[order(pcellsRes7$padj, decreasing = FALSE), ]


gtexLRT <- DESeqDataSet(gtex, design = ~ Plaque + CMC + END + SMC + MAC + ADI) %>% 
  DESeq(test = "LRT", reduced = ~ Plaque + CMC + MAC + SMC + ADI)
lrtRes <- results(gtexLRT, contrast = c("Plaque", 4, 0))
lrtRes$HGNC <- cellResults$HGNC
lrtRes <- lrtRes[,c(7,1:6)]
lrtResSort <- lrtRes[order(lrtRes$padj, decreasing = FALSE),]

toastProp <- t(res_v6$H) %>% data.frame() %>% mutate(Sample = cac$`SRR Number`[specSRR])
rownames(toastProp) <- cac$`SRR Number`[specSRR]

ciberDF <- pivot_longer(ciberOut[,c(1,8,2,3,4,6,7,5)], cols = 2:8, names_to = "Celltype", values_to = "propCIBER")
toastDF <- pivot_longer(toastProp, cols = 1:7, names_to = "Celltype", values_to = "propTOAST") %>% mutate(propCIBER = ciberDF$propCIBER)

symbolsV6 <- vector()
for(i in 1:nrow(rse_coronary)){
  test <- rownames(rse_coronary)[i] %in% gene_tab[,1]
  symbolsV6[i] <- ifelse(test, which(gene_tab[,1] == rownames(rse_coronary)[i]), NA)
}


