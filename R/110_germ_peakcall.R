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

inDir.ArchRProject <- "processed/10_subsetGerm/10_proj_germ"
outDir.ArchRProject <- "processed/20_subsetGerm_peak/20_proj_germ_peak"

outDir.proc <- "processed/20_subsetGerm_peak"
outDir.fig <- "figures/20_subsetGerm_peak"
outDir.res <- "results/20_subsetGerm_peak"
out.table <- "metadata_germ"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

remapClust <- c(
    "C1" = "SSC_1", 
    "C2" = "SSC_2", 
    "C3" = "Prog", 
    "C5" = "Diff", 
    "C6" = "preLep",
    "C8" = "Lep",
    "C9" = "Zyg",
    "C7" = "eP",
    "C4" = "mlP_D",
    "C10" = "eRS", 
    "C11" = "mRS", 
    "C12" = "lRS", 
    "C13" = "ES"
)

labelNew <- mapLabels(names(remapClust), oldLabels = names(remapClust), newLabels = remapClust)
proj$ClustersGerm <- mapLabels(
    proj$Clusters_Germ, 
    newLabels = labelNew, 
    oldLabels = names(remapClust))

#----------------------------------------------------------------------------------------------------
############### Calling peaks ###############
#----------------------------------------------------------------------------------------------------

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersGerm")

pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "ClustersGerm", 
    pathToMacs2 = pathToMacs2)

proj <- addPeakMatrix(proj)

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
