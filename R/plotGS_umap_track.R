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

inDir.ArchRProject <- "processed/02_filtering/01_proj_filtered"

markerGenes <- read.table("gene_sets/markerGenes.txt")[["V1"]]
proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
proj <- addImputeWeights(proj, reducedDims = "IterativeLSI")

#----------------------------------------------------------------------------------------------------
############### Plotting ###############
#----------------------------------------------------------------------------------------------------

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(
    p, 
    name = "Plot-UMAP-RNA-GeneScore", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 8, height = 8)


p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 10000,
    downstream = 10000)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)