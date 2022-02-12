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

inDir.ArchRProject <- "processed/30_subsetGerm_p2g/07_proj_germ_p2g"

outDir.fig <- "figures/plots"
mkdir(outDir.fig)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

markerGenes <- read.table("gene_sets/positiveTFs_germ.txt")[["V1"]]

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]
cluster_name <- "ClustersGerm"

#----------------------------------------------------------------------------------------------------
############### Heatmap of positive regulators ###############
#----------------------------------------------------------------------------------------------------


motif_mat <- getMatrixHeatmap(
    ArchRProj = proj,
    ClustersName = cluster_name,
    MatrixName = "cisBPchromVar_Germ",
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

pdf(file.path(outDir.fig, "Heatmap_PositiveTFs_Germ.pdf"))
Heatmap(motif_mat, cluster_columns = FALSE, cluster_rows = TRUE, paletteContinuous(set = "solarExtra")) + 
Heatmap(gs_mat, cluster_columns = FALSE, cluster_rows = TRUE, paletteContinuous(set = "blueYellow")) +
Heatmap(ge_mat, cluster_columns = FALSE, cluster_rows = TRUE, paletteContinuous(set = "coolwarm"))
dev.off()

# split_cols = factor(cellOrder)
# pdf("Heatmap_PositiveTFs_Germ.pdf")
# Heatmap(
#     motif_mat, 
#     cluster_columns = FALSE, column_split = split_cols, cluster_rows = TRUE, 
#     show_column_names = FALSE, border = "#404040", paletteContinuous(set = "solarExtra")) + 
# Heatmap(gs_mat, 
#     cluster_columns = FALSE, column_split = split_cols, cluster_rows = TRUE, 
#     show_column_names = FALSE, border = "#404040", paletteContinuous(set = "blueYellow")) +
# Heatmap(ge_mat, 
#     cluster_columns = FALSE, column_split = split_cols, cluster_rows = TRUE, 
#     paletteContinuous(set = "coolwarm"))
# dev.off()