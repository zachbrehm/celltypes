makePairs <- function(data) 
{
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
               x = data[, xcol], y = data[, ycol], data)
  }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))
  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
    data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
  }))
  list(all=all, densities=densities)
}

# expand iris data frame for pairs plot
gg1 = makePairs(iris[,-5])
gg2 <- makePairs(data.frame(celltypeProportions))

mega_cell <- data.frame(gg2$all, Plaque = rep(cac$`Plaque Type`, length = nrow(gg2$all)))
mega_cell$Plaque <- factor((mega_cell$Plaque > 0)*1)

# new data frame mega iris
mega_iris = data.frame(gg1$all, Species=rep(iris$Species, length=nrow(gg1$all)))
palettes <- c("#0474BA", "#F17720")
# pairs plot
ggplot(mega_cell, aes_string(x = "x", y = "y")) + 
  facet_grid(xvar ~ yvar, scales = "free") + 
  geom_point(aes(colour=Plaque), na.rm = TRUE, alpha=0.75) + 
  stat_density(aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), 
               data = gg2$densities, position = "identity", 
               colour = "grey20", geom = "line") + 
  scale_color_manual(values = palettes)

