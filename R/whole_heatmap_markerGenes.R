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

inDir.ArchRProject <- "processed/09_trajectory/09_proj_final"

outDir.proc <- "processed/whole"
outDir.fig <- "figures/whole"
outDir.res <- "results/whole"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
markerGenes <- read.table("gene_sets/markerGenes.txt")[["V1"]]
GermTrajectory <- cellOrder[1:10]
cluster_name <- "ClustersWhole"



#----------------------------------------------------------------------------------------------------
############### Heatmap of gene score ###############
#----------------------------------------------------------------------------------------------------


##### Heatmap of gene score to asses quality of cell clusters
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = cluster_name,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
proj <- addImputeWeights(proj, reducedDims = "IterativeLSI")

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE)


pdf(file.path(outDir.fig, "GeneScores-Marker-Heatmap-Whole.pdf"))
ComplexHeatmap::draw(
    heatmapGS, 
    heatmap_legend_side = "bot", 
    annotation_legend_side = "bot",
    row_order = cellOrder)
dev.off()


#----------------------------------------------------------------------------------------------------
############### Heatmap of gene expression ###############
#----------------------------------------------------------------------------------------------------

markersGE <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = cluster_name,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

markerList <- getMarkers(markersGE, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGE <- plotMarkerHeatmap(
    seMarker = markersGE,
    cutOff = "FDR <= 0.01 & Log2FC >= 2",
    pal = ArchRPalettes$coolwarm,
    labelMarkers = markerGenes,
    transpose = TRUE)

pdf(file.path(outDir.fig, "GeneExpression-Marker-Heatmap-Whole.pdf"))
ComplexHeatmap::draw(
    heatmapGE, 
    heatmap_legend_side = "bot", 
    annotation_legend_side = "bot",
    row_order = cellOrder)
dev.off()