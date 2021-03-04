## downloading gtex data with recount3
library(recount3)

## get table of available projects from recount3
human_projects <- available_projects()

## gtex data in recount3 is stored according to tissue type
## we want the blood vessel data in this case
proj_info <- subset(human_projects,
                    project == "BLOOD_VESSEL" & project_type == "data_sources"
                    )

## make the summarized experiment object from the blood vessel gtex data
rse_blood_vessel <- create_rse(proj_info)

## we only want the coronary artery samples, stored in gtex.smtsd as "Artery - Coronary"
rse_coronary <- rse_blood_vessel[,rse_blood_vessel$gtex.smtsd == "Artery - Coronary"]

## transform counts as with the reference data
assay(rse_coronary, "counts") <- transform_counts(rse_coronary)

## save rse for later analysis
saveRDS(rse_coronary, file = "data/rse_coronary.Rds")
