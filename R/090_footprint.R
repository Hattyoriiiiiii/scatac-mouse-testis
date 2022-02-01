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

inDir.ArchRProject <- "processed/07_p2g/07_proj_p2g"

outDir.proc <- "processed/09_footprint"
outDir.fig <- "figures/09_footprint"
outDir.res <- "results/09_footprint"
# out.table <- "metadata_positiveTFs"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersWhole")s

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cluster_name <- "ClustersWhole"

df <- read.csv("results/08_positiveRegulators/metadata_positiveTFs.csv")
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
    plotName = "Footprints-Subtract-Bias-PositiveTFs",
    addDOC = FALSE,
    smoothWindow = 5,
    force = TRUE)
