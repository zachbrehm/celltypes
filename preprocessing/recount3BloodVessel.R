library(recount3)
human_projects <- available_projects()
proj_info <- subset(human_projects,
                    project == "BLOOD_VESSEL" & project_type == "data_sources"
                    )
rse_blood_vessel <- create_rse(proj_info)
rse_coronary <- rse_blood_vessel[,rse_blood_vessel$gtex.smtsd == "Artery - Coronary"]
assay(rse_coronary, "counts") <- transform_counts(rse_coronary)
saveRDS(rse_coronary, file = "data/rse_coronary_recount3.Rds")
