library(DESeq2)
library(tidyverse)
library(Rtsne)

## this script is to perform dimension reduction of gene expression data
## for plotting using the cmd method and tsne method

## read in batch corrected gene expression data
rse <- readRDS(file = "data/ruv_reference_2021_v2.Rds")

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

gg_cmd <- ggplot(cmDf, aes(x = V1, y = V2, text = Project, color = Celltype)) + 
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
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line("grey30", linetype = "dotted"),
        panel.grid.minor = element_line("grey30", linetype = "dotted"),
        legend.key = element_rect(fill = "white"),
        text = element_text(family = "Source Sans Pro"))

ggsave(filename = "preprocessing/gg_cmd.svg", plot = gg_cmd, 
       device = "svg", dpi = 600, width = 11, height = 8.5, units = "in")

## run tsne with the previously created distance matrix
cellTsne <- Rtsne(dist, is_distance = TRUE, perplexity = 40)

## make tsne dataframe and plot
tsneDf <- data.frame(X = cellTsne$Y[,1], Y = cellTsne$Y[,2], Z = rse$celltype)

gg_tsne <- ggplot(tsneDf, aes(x = X, y = Y, color = Z)) + 
  geom_point() + 
  labs(title = "Dimension Reduction with TSNE", 
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
                        begin = 0, end = 0.95) + ## with end = 1 the yellow is a bit too bright for a white background
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line("grey30", linetype = "dotted"),
        panel.grid.minor = element_line("grey30", linetype = "dotted"),
        legend.key = element_rect(fill = "white"),
        text = element_text(family = "Source Sans Pro"))

ggsave(filename = "preprocessing/gg_tsne.svg", plot = gg_tsne, 
       device = "svg", dpi = 600, width = 10.5, height = 8, units = "in")
