suppressPackageStartupMessages({
    library(ArchR)
    library(pheatmap)
    library(tidyverse)
})

#----------------------------------------------------------------------------------------------------
############### Set up ###############
#----------------------------------------------------------------------------------------------------

set.seed(42)
addArchRThreads(threads = 32)
addArchRGenome("Mm10")

inDir.ArchRProject <- "processed/04_Integration_Co/04_proj_CoInteg"
outDir.ArchRProject <- "processed/05_peakcall/05_proj_peakcall"

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

#----------------------------------------------------------------------------------------------------
############### Calling peaks ###############
#----------------------------------------------------------------------------------------------------

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2)

proj <- addPeakMatrix(proj)
saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)