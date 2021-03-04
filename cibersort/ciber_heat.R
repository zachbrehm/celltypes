library(tidyverse)

## read in proportion estimates from cibersort
ciberOut <- read.table(file = "cibersort/output/ciberOut_recount3.txt", header = TRUE)

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
