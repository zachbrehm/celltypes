library(PROPER)

simOpts1 <- RNAseq.SimOptions.2grp(ngenes = 11954, p.DE = 0.01, lOD = "cheung",
                                  lBaselineExpr = "cheung")
simOpts5 <- RNAseq.SimOptions.2grp(ngenes = 11954, p.DE = 0.05, lOD = "cheung",
                                   lBaselineExpr = "cheung")
simOpts10 <- RNAseq.SimOptions.2grp(ngenes = 11954, p.DE = 0.1, lOD = "cheung",
                                   lBaselineExpr = "cheung")

simOptsList <- list(simOpts1, simOpts5, simOpts10)

simResCheungV6.01 <- runSims(Nreps = c(56, 56, 56, 56, 56), Nreps2 = c(19, 35, 11, 10, 84), sim.opts = simOpts1, DEmethod = "DESeq2", nsims = 50)
simResCheungV8.01 <- runSims(Nreps = c(94, 94, 94, 94, 94), Nreps2 = c(32, 51, 18, 28, 146), sim.opts = simOpts1, DEmethod = "DESeq2", nsims = 50)

simResCheungV6.05 <- runSims(Nreps = c(56, 56, 56, 56, 56), Nreps2 = c(19, 35, 11, 10, 84), sim.opts = simOpts5, DEmethod = "DESeq2", nsims = 50)
simResCheungV8.05 <- runSims(Nreps = c(94, 94, 94, 94, 94), Nreps2 = c(32, 51, 18, 28, 146), sim.opts = simOpts5, DEmethod = "DESeq2", nsims = 50)

simResCheungV6.10 <- runSims(Nreps = c(56, 56, 56, 56, 56), Nreps2 = c(19, 35, 11, 10, 84), sim.opts = simOpts10, DEmethod = "DESeq2", nsims = 50)
simResCheungV8.10 <- runSims(Nreps = c(94, 94, 94, 94, 94), Nreps2 = c(32, 51, 18, 28, 146), sim.opts = simOpts10, DEmethod = "DESeq2", nsims = 50)

simV6 <- function(opts){
  runSims(Nreps = c(56, 56, 56, 56, 56), Nreps2 = c(19, 35, 11, 10, 84), sim.opts = opts, DEmethod = "DESeq2", nsims = 50)
}

simV8 <- function(opts){
  runSims(Nreps = c(94, 94, 94, 94, 94), Nreps2 = c(32, 51, 18, 28, 146), sim.opts = opts, DEmethod = "DESeq2", nsims = 50)
}

simResV6 <- mclapply(simOptsList, simV6, mc.cores = 3)
simResV8 <- mclapply(simOptsList, simV8, mc.cores = 3)

powersV6_lfc585 <- lapply(simResV6, comparePower, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 0.585, filter.by = "expr") 
powersV6_lfc1 <- lapply(simResV6, comparePower, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 1) 
powersV8_lfc585 <- lapply(simResV8, comparePower, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 0.585) 
powersV8_lfc1 <- lapply(simResV8, comparePower, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 1) 

powersV6_eff585 <- lapply(simResV6, comparePower, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 0.585) 
powersV6_eff1 <- lapply(simResV6, comparePower, alpha.type = "fdr", alpha.nominal = 0.1,   stratify.by = "dispersion", target.by = "effectsize", delta = 1) 
powersV8_eff585 <- lapply(simResV8, comparePower, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 0.585) 
powersV8_eff1 <- lapply(simResV8, comparePower, alpha.type = "fdr", alpha.nominal = 0.1,   stratify.by = "dispersion", target.by = "effectsize", delta = 1) 

powersList1 <- list(summaryPower(powersV6_lfc585[[1]]), 
                    summaryPower(powersV6_lfc585[[2]]), 
                    summaryPower(powersV6_lfc585[[3]]), 
                    summaryPower(powersV8_lfc585[[1]]), 
                    summaryPower(powersV8_lfc585[[2]]), 
                    summaryPower(powersV8_lfc585[[3]]),
                    summaryPower(powersV6_lfc1[[1]]),
                    summaryPower(powersV6_lfc1[[2]]),
                    summaryPower(powersV6_lfc1[[3]]),
                    summaryPower(powersV8_lfc1[[1]]),
                    summaryPower(powersV8_lfc1[[2]]),
                    summaryPower(powersV8_lfc1[[3]]),
                    summaryPower(powersV6_eff585[[1]]),
                    summaryPower(powersV6_eff585[[2]]),
                    summaryPower(powersV6_eff585[[3]]),
                    summaryPower(powersV8_eff585[[1]]),
                    summaryPower(powersV8_eff585[[2]]),
                    summaryPower(powersV8_eff585[[3]]),
                    summaryPower(powersV6_eff1[[1]]),
                    summaryPower(powersV6_eff1[[2]]),
                    summaryPower(powersV6_eff1[[3]]),
                    summaryPower(powersV8_eff1[[1]]),
                    summaryPower(powersV8_eff1[[2]]),
                    summaryPower(powersV8_eff1[[3]]))



ssPlaque <- rep(c(rep(c("V6 2 vs 0", "V6 4 vs 0", "V6 5 vs 0", "V6 8 vs 0", "V6 All vs 0"), 3), 
              rep(c("V8 2 vs 0", "V8 4 vs 0", "V8 5 vs 0", "V8 8 vs 0", "V8 All vs 0"), 3)),4)

pDE <- rep(c(rep(1, 5), rep(5, 5), rep(10, 5)), 8)

targetDelta <- c(rep("lfc = 0.585", 30), rep("lfc = 1", 30), rep("eff = 0.585", 30), rep("eff = 1", 30))

powerDf <- do.call(rbind, powersList1) %>% data.frame() %>% mutate("Percent DE" = pDE, "Plaque" = ssPlaque, "Target" = targetDelta) %>% select("Percent DE", "Plaque", "Target", everything()) %>% select(-"FDC")
colnames(powerDf) <- c("PercentDE", "Test", "Target", "SS1", "SS2", "Nominal FDR", "Actual FDR", "Marginal Power", "Avg TD", "Avg FD")

powersCheungV6 <- comparePower(simOutput = simResCheungV6, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 0.585)
summaryPower(powersCheungV6)
powersCheungV8 <- comparePower(simOutput = simResCheungV8, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 0.585)
summaryPower(powersCheungV8)

plotAll(powersMAQCV8)

powersCheungV6Eff <- comparePower(simOutput = simResCheungV6, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 1)
summaryPower(powersCheungV6Eff)
powersCheungV8Eff <- comparePower(simOutput = simResCheungV8, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 1)
summaryPower(powersCheungV8Eff)

powersGiladV6 <- comparePower(simOutput = simResGiladV6, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 1)
summaryPower(powersGiladV6)
powersGiladV8 <- comparePower(simOutput = simResGiladV8, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 1)
summaryPower(powersGiladV8)

powersGiladV6Eff <- comparePower(simOutput = simResGiladV6, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 1)
summaryPower(powersGiladV6Eff)
powersGiladV8Eff <- comparePower(simOutput = simResGiladV8, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 1)
summaryPower(powersGiladV8Eff)

powersMAQCV6 <- comparePower(simOutput = simResMAQCV6, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 1)
summaryPower(powersMAQCV6)
powersMAQCV8 <- comparePower(simOutput = simResMAQCV8, alpha.type = "fdr", alpha.nominal = 0.1, target.by = "lfc", delta = 1)
summaryPower(powersMAQCV8)

powersMAQCV6Eff <- comparePower(simOutput = simResMAQCV6, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 1)
summaryPower(powersMAQCV6Eff)
powersMAQCV8Eff <- comparePower(simOutput = simResMAQCV8, alpha.type = "fdr", alpha.nominal = 0.1, stratify.by = "dispersion", target.by = "effectsize", delta = 1)
summaryPower(powersMAQCV8Eff)

powersList1 <- list(summaryPower(powersCheungV6), 
                   summaryPower(powersCheungV8), 
                   summaryPower(powersGiladV6), 
                   summaryPower(powersGiladV8), 
                   summaryPower(powersMAQCV6), 
                   summaryPower(powersMAQCV8))

refDat <- c(rep("Cheung", 4),
            rep("Gilad", 4),
            rep("MAQC", 4))

ssPlaque <- rep(c("V6 4 vs 0", "V6 All vs 0", "V8 4 vs 0", "V8 All vs 0"), 3)

powerDf1 <- do.call(rbind, powersList1) %>% data.frame() %>% mutate("Reference" = refDat, "Plaque" = ssPlaque) %>% select("Reference", "Plaque", everything())
colnames(powerDf1) <- c("Reference", "Plaque comparison for SS", "# Plaque Group 1", "# Plaque Group 2", "Nominal FDR", "Actual FDR", "Marginal Power", "Avg # of TD", "Avg # of FD", "FDC")
