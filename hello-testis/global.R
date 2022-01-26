library(ArchR)
library(tidyverse)

set.seed(42)

print("Please wait...")
proj <- loadArchRProject(path = "/work/Save-proj_traj/", showLogo = FALSE)
proj <- addImputeWeights(proj)

motifPositions <- getPositions(proj, name = "Motif_cisbp")
genes <- proj@geneAnnotation$genes$symbol %>% unique()
GermTrajectory <- c("Undiff_1", "Undiff_2", "Diff", "Lep", "Zyg_eP", "Pac_Dip", "RS_1", "RS_2", "RS_3", "ES")
cellOrder <- c("Undiff_1", "Undiff_2", "Diff", "Lep", "Zyg_eP", "Pac_Dip", "RS_1", "RS_2", "RS_3", "ES", "Sertoli", "Soma")
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersAnno", force=TRUE)


print(getAvailableMatrices(proj))
