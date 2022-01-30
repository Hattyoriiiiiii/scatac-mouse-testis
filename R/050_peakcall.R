suppressPackageStartupMessages({
    library(ArchR)
    library(pheatmap)
    library(tidyverse)
})

#----------------------------------------------------------------------------------------------------
############### Import functions ###############
#----------------------------------------------------------------------------------------------------

source("R/functions.R")

#----------------------------------------------------------------------------------------------------
############### Set up ###############
#----------------------------------------------------------------------------------------------------

set.seed(42)
addArchRThreads(threads = 32)
addArchRGenome("Mm10")

inDir.ArchRProject <- "processed/04_Integration_Co/04_proj_CoInteg"
outDir.ArchRProject <- "processed/05_peakcall/05_proj_peakcall"

outDir.proc <- "processed/05_peakcall"
mkdir(outDir.proc)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

remapClust <- c(
    "C1" = "DM", 
    "C2" = "mlP", 
    "C3" = "preLep", 
    "C4" = "ZeP",
    "C5" = "Lep", 
    "C6" = "Sertoli",
    "C7" = "Soma",
    "C8" = "SSC",
    "C9" = "Diff",
    "C10" = "Prog", 
    "C11" = "eRS", 
    "C12" = "mRS", 
    "C13" = "lRS", 
    "C14" = "ES" 
)

labelNew <- mapLabels(names(remapClust), oldLabels = names(remapClust), newLabels = remapClust)
proj$ClustersWhole <- mapLabels(
    proj$Clusters, 
    newLabels = labelNew, 
    oldLabels = names(remapClust))

#----------------------------------------------------------------------------------------------------
############### Calling peaks ###############
#----------------------------------------------------------------------------------------------------

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersWhole")

pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "ClustersWhole", 
    pathToMacs2 = pathToMacs2)

proj <- addPeakMatrix(proj)
saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)