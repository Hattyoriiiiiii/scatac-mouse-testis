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
outDir.ArchRProject <- "processed/10_subsetGerm/10_proj_germ"

outDir.proc <- "processed/10_subsetGerm"
outDir.fig <- "figures/10_subsetGerm"
outDir.res <- "results/10_subsetGerm"
out.table <- "metadata_germ"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cellOrder_germ <- cellOrder[1:10]

markerGenes <- read.table("gene_sets/RepresentativeGerm.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

proj_germ <- subsetCells(proj, cellNames = proj[proj$ClustersWhole %in% cellOrder_germ]$cellNames)

##### Iterative Latent Semantic Indexing -----
proj_germ <- addIterativeLSI(
    ArchRProj = proj_germ,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI_Germ", 
    iterations = 5, 
    clusterParams = list( 
        resolution = c(0.1, 0.2, 0.4, 0.8), 
        sampleCells = 10000, 
        n.start = 10), 
    varFeatures = 30000, 
    dimsToUse = 1:50,
    force = TRUE)


# Clustering using Seurat’s FindClusters() function -----
proj_germ <- addClusters(
    input = proj_germ,
    reducedDims = "IterativeLSI_Germ",
    method = "Seurat",
    name = "Clusters_Germ",
    resolution = 1.6,
    force = TRUE)

# UMAP embedding -----
proj_germ <- addUMAP(
    ArchRProj = proj_germ, 
    reducedDims = "IterativeLSI_Germ", 
    name = "UMAP_Germ", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE)

proj_germ <- addImputeWeights(proj_germ, reducedDims = "IterativeLSI_Germ")

#----------------------------------------------------------------------------------------------------
###############  ###############
#----------------------------------------------------------------------------------------------------

pdf(file.path(outDir.fig, "Germ-UMAP.pdf"))
qc_plot(
    proj_germ, 
    qc_metrix = c("Sample", "Clusters_Germ", "log10(nFrags)", "DoubletEnrichment", "TSSEnrichment"), 
    embedding = "UMAP_Germ")
dev.off()

proj_germ@cellColData %>%
    as.data.frame() %>%
    write.table(
        file = paste0(file.path(outDir.res, out.table), ".csv"),
        quote = FALSE,
        sep = ",",
        row.names = TRUE,
        col.names = TRUE)


##### Gene score -----
markersGS <- getMarkerFeatures(
    ArchRProj = proj_germ, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters_Germ",
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
    ArchRProj = proj_germ, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = "Clusters_Germ",
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
    ArchRProj = proj_germ, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)