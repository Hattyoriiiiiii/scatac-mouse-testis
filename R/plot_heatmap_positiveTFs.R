suppressPackageStartupMessages({
    library(ArchR)
    library(pheatmap)
    library(tidyverse)
    library(ComplexHeatmap)
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

outDir.fig <- "figures/plots"
mkdir(outDir.fig)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

markerGenes <- read.table("gene_sets/positiveTFs_whole.txt")[["V1"]]

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cluster_name <- "ClustersWhole"

#----------------------------------------------------------------------------------------------------
############### Heatmap of positive regulators ###############
#----------------------------------------------------------------------------------------------------


motif_mat <- getMatrixHeatmap(
    ArchRProj = proj,
    ClustersName = cluster_name,
    MatrixName = "cisBPchromVar",
    markerGenes = markerGenes,
    cellOrder = cellOrder,
    motif = TRUE
)

gs_mat <- getMatrixHeatmap(
    ArchRProj = proj,
    ClustersName = cluster_name,
    MatrixName = "GeneScoreMatrix",
    markerGenes = markerGenes,
    cellOrder = cellOrder,
    motif = FALSE
)

ge_mat <- getMatrixHeatmap(
    ArchRProj = proj,
    ClustersName = cluster_name,
    MatrixName = "GeneIntegrationMatrix",
    markerGenes = markerGenes,
    cellOrder = cellOrder,
    motif = FALSE
)

pdf(file.path(outDir.fig, "Heatmap_PositiveTFs.pdf"))
Heatmap(motif_mat, cluster_columns = FALSE, cluster_rows = TRUE, paletteContinuous(set = "solarExtra")) + 
Heatmap(gs_mat, cluster_columns = FALSE, cluster_rows = TRUE, paletteContinuous(set = "blueYellow")) +
Heatmap(ge_mat, cluster_columns = FALSE, cluster_rows = TRUE, paletteContinuous(set = "coolwarm"))
dev.off()