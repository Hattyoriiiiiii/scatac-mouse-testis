library(ArchR)
# library(tidyverse)

set.seed(42)

print("Please wait...")
proj <- loadArchRProject(path = "/work/Save-proj_traj/", showLogo = FALSE)
# proj <- addImputeWeights(proj)

motifPositions <- readRDS("/work/hello-testis/motifPositions.rds")
genes <- readRDS("/work/hello-testis/genes.rds")
GermTrajectory <- c("Undiff_1", "Undiff_2", "Diff", "Lep", "Zyg_eP", "Pac_Dip", "RS_1", "RS_2", "RS_3", "ES")
cellOrder <- c("Undiff_1", "Undiff_2", "Diff", "Lep", "Zyg_eP", "Pac_Dip", "RS_1", "RS_2", "RS_3", "ES", "Sertoli", "Soma")
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersAnno", force=TRUE)


print(getAvailableMatrices(proj))
