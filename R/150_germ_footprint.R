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

inDir.ArchRProject <- "processed/30_subsetGerm_p2g/07_proj_germ_p2g"

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersWhole")s

cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]
cluster_name <- "ClustersGerm"

df <- read.csv("results/14_positiveRegulators/metadata_positiveTFs.csv")
motifs <- df[df$TFRegulator == "YES",]$GeneIntegrationMatrix_name

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

motifPositions <- getPositions(proj)
markerMotifs <- unlist(lapply(paste0("^", motifs, "_"), function(x) grep(x, names(motifPositions), value = TRUE)))

### TF footprinting -----
seFoot <- getFootprints(
    ArchRProj = proj, 
    positions = motifPositions[markerMotifs], 
    groupBy = cluster_name)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias-PositiveTFs-Germ",
    addDOC = FALSE,
    smoothWindow = 5,
    force = TRUE)
