## Make a reference dataset for CIBERSORT with recount3
library(recount3)

## have an object containing the projects that you want to download
## we are getting data from the SRA database
## have project names in the SRP naming convention, i.e., SRP009867
projects <- readRDS("data/projects.Rds")

## get table of available projects for download
human_projects <- available_projects()

## subset initial project names to those that are available in recount3
projects <- projects[projects %in% human_projects$project]

## select first project to initialize rse object
project <- subset(human_projects, project == projects[1] & project_type == "data_sources")

## download first project and store as rse
rse <- create_rse(project_info = project)

## remove first project from projects vector
projects <- projects[-1]

## loop over remaining projects
## first isolate one project for download since recount3 takes one project at a time
## create temporary rse object
## cbind rse and tmp_rse
for(i in 1:length(projects)){
  project <- subset(human_projects, project == projects[i] & project_type == "data_sources")
  tmp_rse <- create_rse(project_info = project)
  rse <- cbind(rse, tmp_rse)
}

## Have some object available that identifies the celltypes of each sample in the reference data
celltypes <- readRDS("data/celltypes.Rds")

## assign celltypes to samples
rse$celltype <- celltypes

## recount3 gives the raw count data, we need to transform these counts to scale according to a genes abundance
assay(rse, "counts") <- transform_counts(rse)

## save rse reference object
saveRDS(rse, file = "data/rse_reference.Rds")
