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

inDir.ArchRProject <- "processed/10_subsetGerm/10_proj_germ"

markerGenes <- read.table("gene_sets/markerGenes.txt")[["V1"]]
proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
proj <- addImputeWeights(proj, reducedDims = "IterativeLSI_Germ")

#----------------------------------------------------------------------------------------------------
############### Plotting ###############
#----------------------------------------------------------------------------------------------------

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAP_Germ",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(
    p, 
    name = "Plot-UMAP-RNA-GeneScore", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 8, height = 8)


p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters_Germ", 
    geneSymbol = markerGenes, 
    upstream = 10000,
    downstream = 10000)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-Representative.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)