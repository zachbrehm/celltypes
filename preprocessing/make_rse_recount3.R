projects <- unique(rse_6_2_20$project)
projects <- projects[-(which(projects %in% c("Halushka", "DRP002851")))]

library(recount3)

human_projects <- available_projects()

projects <- projects[projects %in% human_projects$project]

project <- subset(human_projects, project == projects[1] & project_type == "data_sources")

rse <- create_rse(project_info = project)
projects <- projects[-1]
for(i in 1:length(projects)){
  project <- subset(human_projects, project == projects[i] & project_type == "data_sources")
  tmp_rse <- create_rse(project_info = project)
  rse <- cbind(rse, tmp_rse)
}

samples <- which(colnames(rse) %in% colnames(rse_6_2_20))

rse <- rse[,samples]
columnsSums <- colSums(assay(rse))
rse <- rse[,-which(columnsSums < 10)]
rse_old <- rse_6_2_20[,colnames(rse_6_2_20) %in% colnames(rse)]
map <- match(colnames(rse_old), colnames(rse))
rse <- rse[,map]

rse$celltype <- rse_old$celltype

assay(rse, "counts") <- transform_counts(rse)

saveRDS(rse, file = "data/rse_reference_2021_v1.Rds")
