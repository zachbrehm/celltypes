#load(file = "rseData_01_06_2020.RData")
#rse <- rseData_01_06_2020
### look just in DRO001797
#ind <- which(colData(rse)$project == "DRP001797")
#
#runs <- rse$run[ind]
#
#samples <- rse$sample[ind] %>% unique()
#
### matrix to contain rowSums across replicates
#sumMat <- matrix(nrow = nrow(rse), ncol = length(samples))
#
### rowSum across columns for replicates from each sample
#for(i in 1:length(samples)){
#  sampInd <- which(rse$sample == samples[i])
#  sumMat[,i] <- rowSums(assay(rse)[,sampInd])
#}
#
#rmVec <- vector()
#keepVec <- vector()
#
#for(i in 1:length(samples)){
#  tmp <- which(rse$sample == samples[i]) %>% sort()
#  rmVec <- c(rmVec, tmp[-1])
#  keepVec <- c(keepVec, tmp[1])
#}
#
#for(i in 1:length(samples)){
#  assay(rse)[,keepVec[i]] <- sumMat[,i]
#}
#
#rse$project[which(is.na(rse$project))] <- "Halushka"
#
#rseData_02_25_2020 <- rse[,-rmVec]
#
#rm(rse)
#
#saveRDS(rseData_02_25_2020, file = "rseData_02_25_2020.Rds")
#
#rseTmp <- readRDS(file = "rse_noFoamFibro.Rds")
#
#removeRuns <- c("SRR307899", "SRR307900", "SRR315297", "SRR315298", "SRR317059", "SRR317060",
#                "SRR317061",
#                "SRR545687",
#                "SRR545688",
#                "SRR1550981",
#                "SRR1550982",
#                "SRR1550983",
#                "SRR1550984",
#                "SRR1550985",
#                "SRR1551019",
#                "SRR1551025",
#                "SRR1551032",
#                "SRR1551039",
#                "SRR1551046",
#                "SRR1551053",
#                "SRR1551060",
#                "SRR1551067",
#                "SRR1551074",
#                "SRR1551081",
#                "SRR1551087",
#                "SRR1551094",
#                "SRR1551107",
#                "SRR1551114",
#                "SRR1609989",
#                "SRR1609990",
#                "SRR1852860",
#                "SRR1852861",
#                "SRR2044053",
#                "SRR2044054",
#                "SRR2044056",
#                "SRR2044057",
#                "SRR2044059",
#                "SRR2044060",
#                "SRR2044062",
#                "SRR2044063",
#                "SRR2044065",
#                "SRR2044066",
#                "SRR2044068",
#                "SRR2044069",
#                "SRR2044071",
#                "SRR2044072",
#                "SRR2044074",
#                "SRR2044075",
#                "SRR2044083",
#                "SRR2044084",
#                "SRR2044086",
#                "SRR2044087",
#                "SRR2044095",
#                "SRR2044096",
#                "SRR2044098",
#                "SRR2044099",
#                "SRR2050392",
#                "SRR2050393",
#                "SRR2050395",
#                "SRR2050396",
#                "SRR2050397",
#                "SRR2050398",
#                "SRR2050400",
#                "SRR2050401",
#                "SRR2050403",
#                "SRR2050405",
#                "SRR2050406",
#                "SRR2050408",
#                "SRR2050409",
#                "SRR2050411",
#                "SRR2050412",
#                "SRR2050414",
#                "SRR2050415",
#                "SRR2050417",
#                "SRR2050418",
#                "SRR2050420",
#                "SRR2050421",
#                "SRR2050423",
#                "SRR2050424",
#                "SRR2050426",
#                "SRR2050427",
#                "SRR2050429",
#                "SRR2050430",
#                "SRR2050431",
#                "SRR2050432",
#                "SRR2050434",
#                "SRR2050435",
#                "SRR2050437",
#                "SRR2050438",
#                "SRR2050440",
#                "SRR2050441",
#                "SRR2050443",
#                "SRR2050444",
#                "SRR2050446",
#                "SRR2050447",
#                "SRR2050448",
#                "SRR2050449",
#                "SRR2050451",
#                "SRR2050452",
#                "SRR2050454",
#                "SRR2050455",
#                "SRR2050457",
#                "SRR2050458",
#                "SRR2050460",
#                "SRR2050461",
#                "SRR2050463",
#                "SRR2050464",
#                "SRR2138073",
#                "SRR2138074",
#                "SRR2138075",
#                "SRR2138076",
#                "SRR2138077",
#                "SRR2138078",
#                "SRR2138079",
#                "SRR2138080",
#                "SRR2138081",
#                "SRR2138082",
#                "SRR2138083",
#                "SRR2138084",
#                "SRR2138085",
#                "SRR2138086",
#                "SRR2138087",
#                "SRR2138088",
#                "SRR2138089",
#                "SRR2138090",
#                "SRR2138091",
#                "SRR2138092",
#                "SRR2138093",
#                "SRR2138094",
#                "SRR2138095",
#                "SRR2138096")
#
#removeInd <- vector(length = length(removeRuns))
#
#for(i in 1:length(removeInd)){
#  removeInd[i] <- which(rseTmp$run == removeRuns[i])
#}
#
#rseTmp <- rseTmp[,-removeInd]
#
#moreRemoveRuns <- c("SRR834938", "SRR834939", "SRR834940", "SRR834941", "SRR834942", "SRR834943", "SRR834944", "SRR834945", "SRR834946", "SRR834947", "SRR834948", "SRR834949", "SRR834950", "SRR834951", "SRR834952", "SRR834953",
#                    "SRR834954", "SRR834955", "SRR834956", "SRR834957", "SRR834958", "SRR834959", "SRR834960", "SRR834961", "SRR834962", "SRR834963", "SRR834964", "SRR834965", "SRR834966", "SRR834967", "SRR834968", "SRR834969",
#                    "SRR834970", "SRR834971", "SRR834972")
#
#moreRemoveInd <- vector(length = length(moreRemoveRuns))
#
#for(i in 1:length(moreRemoveInd)){
#  moreRemoveInd[i] <- which(rseTmp$run == moreRemoveRuns[i])
#}
#
#rseTmp <- rseTmp[,-moreRemoveInd]
#
#lymph_to_macro_runs <- c("SRR1550986",
#                         "SRR1550987",
#                         "SRR1550992",
#                         "SRR1550993",
#                         "SRR1550997",
#                         "SRR1550998",
#                         "SRR1551002",
#                         "SRR1551003",
#                         "SRR1551008",
#                         "SRR1551009",
#                         "SRR1551013",
#                         "SRR1551014",
#                         "SRR1551020",
#                         "SRR1551021",
#                         "SRR1551026",
#                         "SRR1551027",
#                         "SRR1551033",
#                         "SRR1551034",
#                         "SRR1551040",
#                         "SRR1551041",
#                         "SRR1551047",
#                         "SRR1551048",
#                         "SRR1551054",
#                         "SRR1551055",
#                         "SRR1551061",
#                         "SRR1551062",
#                         "SRR1551068",
#                         "SRR1551069",
#                         "SRR1551075",
#                         "SRR1551076",
#                         "SRR1551082",
#                         "SRR1551083",
#                         "SRR1551088",
#                         "SRR1551089",
#                         "SRR1551095",
#                         "SRR1551096",
#                         "SRR1551101",
#                         "SRR1551102",
#                         "SRR1551108",
#                         "SRR1551109",
#                         "SRR1613933",
#                         "SRR1613937",
#                         "SRR1740038",
#                         "SRR1740039",
#                         "SRR1740040",
#                         "SRR1740041",
#                         "SRR1740042",
#                         "SRR1740043",
#                         "SRR1740044",
#                         "SRR1740045",
#                         "SRR1740046",
#                         "SRR1740047",
#                         "SRR1740048",
#                         "SRR1740049",
#                         "SRR1740066",
#                         "SRR1740067",
#                         "SRR1740068",
#                         "SRR1740069",
#                         "SRR1740070",
#                         "SRR1740071",
#                         "SRR1740072",
#                         "SRR1740073",
#                         "SRR1740074",
#                         "SRR1740075",
#                         "SRR1740076",
#                         "SRR1740077")
#
#for(i in 1:length(lymph_to_macro_runs)){
#  tmp <- which(rseTmp$run == lymph_to_macro_runs[i])
#  rseTmp$celltype[tmp] <- "Macrophage"
#}
#
#evenMoreRemove <- c("SRR2910666",
#                    "SRR2910667",
#                    "SRR2939137",
#                    "SRR2939138",
#                    "SRR2939139",
#                    "SRR2939140")
#
#evenMoreInd <- vector(length = length(evenMoreRemove))
#
#for(i in 1:length(evenMoreInd)){
#  evenMoreInd[i] <- which(rseTmp$run == evenMoreRemove[i])
#}
#
#rseTmp <- rseTmp[,-evenMoreInd]
#
#lastRemove <- c("SRR1947490",
#                "SRR1947491",
#                "SRR1947492",
#                "SRR1947494",
#                "SRR1947495",
#                "SRR1947496",
#                "SRR1947497",
#                "SRR1947498",
#                "SRR1947499",
#                "SRR1947500",
#                "SRR1947501",
#                "SRR1947502",
#                "SRR1947503",
#                "SRR1947504",
#                "SRR1947505",
#                "SRR1947506",
#                "SRR1947508",
#                "SRR1947509",
#                "SRR1947510",
#                "SRR1947511",
#                "SRR1947513",
#                "SRR1947514",
#                "SRR1947515",
#                "SRR1947516",
#                "SRR1947517",
#                "SRR1947518",
#                "SRR1947519",
#                "SRR1947520",
#                "SRR1947521",
#                "SRR1947523",
#                "SRR1947524",
#                "SRR1947525",
#                "SRR1947526",
#                "SRR1947527",
#                "SRR1947528",
#                "SRR1947529",
#                "SRR1947530",
#                "SRR1947531",
#                "SRR1947532",
#                "SRR1947533",
#                "SRR1947534",
#                "SRR1947535",
#                "SRR1947536",
#                "SRR1947537",
#                "SRR1947538",
#                "SRR1947539",
#                "SRR1947540",
#                "SRR1947541",
#                "SRR1947542",
#                "SRR1947546",
#                "SRR1947547",
#                "SRR1947548",
#                "SRR1947549",
#                "SRR1947550",
#                "SRR1947551",
#                "SRR1947552",
#                "SRR1947553",
#                "SRR1947554",
#                "SRR1947555",
#                "SRR1947558",
#                "SRR1947559",
#                "SRR1947560",
#                "SRR1947561",
#                "SRR1947562",
#                "SRR1947563",
#                "SRR1947564",
#                "SRR1947565",
#                "SRR1947566",
#                "SRR1947567",
#                "SRR1947568",
#                "SRR1947571",
#                "SRR1947572",
#                "SRR1947573",
#                "SRR1947574",
#                "SRR1947575",
#                "SRR1947576",
#                "SRR1947577",
#                "SRR1947578",
#                "SRR1947579",
#                "SRR1947580",
#                "SRR1947582",
#                "SRR1947583",
#                "SRR1947585",
#                "SRR1947586",
#                "SRR1947587",
#                "SRR1947588",
#                "SRR1947591",
#                "SRR1947592",
#                "SRR1947594",
#                "SRR1947595",
#                "SRR1947597",
#                "SRR1947599",
#                "SRR1947600",
#                "SRR1947601",
#                "SRR1947602",
#                "SRR1947603",
#                "SRR1947604",
#                "SRR1947606",
#                "SRR1947607",
#                "SRR1947608",
#                "SRR1947609",
#                "SRR1947610",
#                "SRR1947611",
#                "SRR1947612",
#                "SRR1947613",
#                "SRR1947614",
#                "SRR1947616",
#                "SRR1947617",
#                "SRR1947618",
#                "SRR1947620",
#                "SRR1947621",
#                "SRR1947622",
#                "SRR1947623",
#                "SRR1947624",
#                "SRR1947625",
#                "SRR1947626",
#                "SRR1947627",
#                "SRR1947628",
#                "SRR1947629",
#                "SRR1947630",
#                "SRR1947631",
#                "SRR1947632",
#                "SRR1947633",
#                "SRR1947634",
#                "SRR1947635",
#                "SRR1947636",
#                "SRR1947637",
#                "SRR1947638",
#                "SRR1947639",
#                "SRR1947640",
#                "SRR1947641",
#                "SRR1947642",
#                "SRR1947643",
#                "SRR1947644")
#
#lastRemoveInd <- vector(length = length(lastRemove))
#
#for(i in 1:length(lastRemoveInd)){
#  lastRemoveInd[i] <- which(rseTmp$run == lastRemove[i])
#}
#
#rseTmp <- rseTmp[,-lastRemoveInd]
#
#macro_to_lymph_runs <- c("SRR534311",
#                         "SRR534312",
#                         "SRR534313",
#                         "SRR534314",
#                         "SRR534315",
#                         "SRR534316",
#                         "SRR534319",
#                         "SRR534320",
#                         "SRR534321",
#                         "SRR534322",
#                         "SRR534323",
#                         "SRR534324")
#
#macro_to_lypmh_ind <- vector(length = length(macro_to_lymph_runs))
#
#for(i in 1:length(macro_to_lymph_runs)){
#  tmp <- which(rseTmp$run == macro_to_lymph_runs[i])
#  rseTmp$celltype[tmp] <- "Lymphocytes"
#}
#
#saveRDS(rseTmp, file = "rse_6_2_20.Rds")

rse <- readRDS(file = "rse_6_2_20.Rds")
adiInd <- which(rse$celltype == "FatCell")
rse <- rse[,-adiInd]

load(file = "rse_gene_adipose_tissue.Rdata")
rse_gene <- scale_counts(rse_gene)
rse_gene <- rse_gene[, which(rse_gene$smtsd == "Adipose - Visceral (Omentum)")]
colData(rse_gene) <- colData(rse_gene)[,1:19]
rse_gene$celltype <- rep("Adipose", ncol(rse_gene))
rse <- cbind(rse, rse_gene)

saveRDS(rse, file = "rse_adipose.Rds")
