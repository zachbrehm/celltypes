library(tidyverse)
library(sva)
library(bladderbatch)
data(bladderdata)
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)
mod <- model.matrix(~as.factor(cancer), data = pheno)
mod0 <- model.matrix(~1, data = pheno)
n.sv <- num.sv(edata, mod, method = "leek")
svobj <- sva(edata, mod, mod0, n.sv)
modSv <- cbind(mod, svobj$sv)
mod0Sv <- cbind(mod0, svobj$sv)
fit <- lmFit(edata, modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),
                         "C2"=c(0,-1,1,rep(0,svobj$n.sv)),
                         "C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts <- contrasts.fit(fit, contrast.matrix)
eb <- eBayes(fitContrasts)
topTableF(eb, adjust = "BH")

batch <- pheno$batch
modcombat <- model.matrix(~1, data = pheno)
combat_edata <- ComBat(dat = edata, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
