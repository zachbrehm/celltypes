library(DESeq2)
library(tidyverse)
library(plotly)
rse <- readRDS(file = "data/ruv_reference_2021_v1.Rds")
dds <- DESeqDataSet(rse, design = ~1)
vst <- varianceStabilizingTransformation(dds)
dist <- dist(t(assay(vst)))
cmd <- cmdscale(dist)
cmDf <- as.data.frame(cmd) %>% 
  rownames_to_column(var = "Sample") %>% 
  mutate(Celltype = rse$celltype, 
         Project = rse$study) %>% 
  select(Sample, 
         Celltype,
         Project,
         everything())
gg <- ggplot(cmDf, aes(x = V1, y = V2, text = Project, color = Celltype)) + 
  geom_point() + labs(title = "Cell Type Clustering", x = "First Dimension", y = "Second Dimension", color = "Cell Type") + 
  theme() + scale_color_viridis_d(labels = c("Cardiac Muscle", "Endothelial", "Erythrocyte", "Adipocyte", "Lymphocyte", "Macrophage", "Smooth Muscle"))
ggplotly(gg)
cellTsne <- Rtsne(dist, is_distance = TRUE, perplexity = 10)
ggplot(data.frame(X = cellTsne$Y[,1], 
                  Y = cellTsne$Y[,2], 
                  Z = rse$celltype), 
       aes(x = X, y = Y, color = Z)) + geom_point() + scale_color_pomological() + theme_pomological()

cellUmap <- umap(d = as.matrix(dist), input = "dist")
ggplot(data.frame(X = cellUmap$layout[,1],
                  Y = cellUmap$layout[,2],
                  Z = rse$celltype),
       aes(x = X, y = Y, color = Z)) + geom_point() + scale_color_pomological()
