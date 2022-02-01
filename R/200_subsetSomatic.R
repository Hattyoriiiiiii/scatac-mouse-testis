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
outDir.ArchRProject <- "processed/20_subsetSomatic/20_proj_somatic"

outDir.proc <- "processed/20_subsetSomatic"
outDir.fig <- "figures/20_subsetSomatic"
outDir.res <- "results/20_subsetSomatic"
out.table <- "metadata_somatic"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cellOrder_somatic <- cellOrder[13:14]

markerGenes <- read.table("gene_sets/markerGenes_somatic.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

proj_soma <- subsetCells(proj, cellNames = proj[proj$ClustersWhole %in% cellOrder_somatic]$cellNames)

##### Iterative Latent Semantic Indexing -----
proj_soma <- addIterativeLSI(
    ArchRProj = proj_soma,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI_Soma", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1, 0.2, 0.4, 0.8), 
        sampleCells = NULL, 
        n.start = 10), 
    varFeatures = 30000, 
    dimsToUse = 1:50,
    force = TRUE)


# Clustering using Seurat’s FindClusters() function -----
proj_soma <- addClusters(
    input = proj_soma,
    reducedDims = "IterativeLSI_Soma",
    method = "Seurat",
    name = "Clusters_Soma",
    resolution = 0.8,
    force = TRUE)

# UMAP embedding -----
proj_soma <- addUMAP(
    ArchRProj = proj_soma, 
    reducedDims = "IterativeLSI_Soma", 
    name = "UMAP_Soma", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE)

proj_soma <- addImputeWeights(proj_soma, reducedDims = "IterativeLSI_Soma")

#----------------------------------------------------------------------------------------------------
###############  ###############
#----------------------------------------------------------------------------------------------------

pdf(file.path(outDir.fig, "Somatic-UMAP.pdf"))
qc_plot(
    proj_soma, 
    qc_metrix = c("Sample", "Clusters_Soma", "log10(nFrags)", "DoubletEnrichment", "TSSEnrichment"), 
    embedding = "UMAP_Soma")
dev.off()

proj_soma@cellColData %>%
    as.data.frame() %>%
    write.table(
        file = paste0(file.path(outDir.res, out.table), ".csv"),
        quote = FALSE,
        sep = ",",
        row.names = TRUE,
        col.names = TRUE)


##### Gene score -----
markersGS <- getMarkerFeatures(
    ArchRProj = proj_soma, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters_Soma",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE)

pdf(file.path(outDir.fig, "GeneScore-Marker-Heatmap.pdf"))
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


##### Gene expression -----
markersGS <- getMarkerFeatures(
    ArchRProj = proj_soma, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = "Clusters_Soma",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 2",
    pal = ArchRPalettes$coolwarm,
    labelMarkers = markerGenes,
    transpose = TRUE)

pdf(file.path(outDir.fig, "GeneExpression-Marker-Heatmap.pdf"))
ComplexHeatmap::draw(
    heatmapGS, 
    heatmap_legend_side = "bot", 
    annotation_legend_side = "bot")
dev.off()


saveArchRProject(
    ArchRProj = proj_soma, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)