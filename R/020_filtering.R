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

inDir.ArchRProject <- "processed/01_qc/01_proj_qc"

outDir.ArchRProject <- "processed/02_filtering/01_proj_filtered"
outDir.proc <- "processed/02_filtering"
outDir.fig <- "figures/02_filtering"
outDir.res <- "results/02_filtering"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

qc.table <- "QC_metadata_after_filtering"
remove.table <- "QC_metadata_removeClusters"

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

markerGenes <- read.table("gene_sets/markerGenes.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Filtering putative doublets / low quality cells ###############
#----------------------------------------------------------------------------------------------------

clustersRemove <- getCellColData(proj) %>%
    as.data.frame() %>% 
    group_by(Clusters_qc) %>% 
    select(TSSEnrichment, nFrags, DoubletEnrichment) %>% 
    summarize(
        med_Doub = median(DoubletEnrichment),
        med_nFrag = median(nFrags),
        Doublet = med_Doub > 1
    ) %>% 
    filter(Doublet == TRUE) %>% 
    select(Clusters_qc) %>% 
    as.vector

getCellColData(proj) %>%
    as.data.frame() %>% 
    group_by(Clusters_qc) %>% 
    select(TSSEnrichment, nFrags, DoubletEnrichment) %>% 
    summarize(
        med_Doub = median(DoubletEnrichment),
        med_nFrag = median(nFrags),
        Doublet = med_Doub > 1) %>% 
    write.table(
        file = paste0(file.path(outDir.res, remove.table), ".csv"),
        quote = FALSE,
        sep = ",",
        row.names = TRUE,
        col.names = TRUE)

proj_sub <- subsetCells(
    proj, 
    cellNames = proj[proj$Clusters_qc %ni% clustersRemove$Clusters_qc &
                       proj$nFrags < 25000 & 
                       proj$DoubletEnrichment < 2.5]$cellNames)


#----------------------------------------------------------------------------------------------------
############### Dimensionality Reduction, Clustering ###############
#----------------------------------------------------------------------------------------------------

##### Iterative Latent Semantic Indexing -----
proj_sub <- addIterativeLSI(
    ArchRProj = proj_sub,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1, 0.2, 0.4, 0.8), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 30000, 
    dimsToUse = 1:50,
    excludeChr = c("chrX", "chrY", "chrMT"),
    force = TRUE)

# Clustering using Seurat’s FindClusters() function -----
proj_sub <- addClusters(
    input = proj_sub,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 1.2,
    force = TRUE)

# UMAP embedding -----
proj_sub <- addUMAP(
    ArchRProj = proj_sub, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE)


#----------------------------------------------------------------------------------------------------
############### Quality Assessment ###############
#----------------------------------------------------------------------------------------------------

pdf(file.path(outDir.fig, "QC-subsetting.pdf"))
qc_plot(
    proj_sub, 
    qc_metrix = c("Sample", "Clusters", "log10(nFrags)", "DoubletEnrichment", "TSSEnrichment"),
    embedding = "UMAP")
dev.off()


pdf(file.path(outDir.fig, "QC-after-subsetting.pdf"))
p <- ggPoint(
        x=proj_sub$nFrags,
        y=proj_sub$DoubletEnrichment,
        colorDensity = T,
        rastr = T) +
    geom_hline(yintercept = 3, col="red", lty="dashed") +
    geom_vline(xintercept = 30000, col="red", lty="dashed") +
    xlab("nFragments") +
    ylab("DoubletEnrichment")
print(p)
dev.off()

# save the quality metrix
getCellColData(proj_sub) %>%
    as.data.frame() %>% 
    group_by(Clusters) %>% 
    select(TSSEnrichment, nFrags, DoubletEnrichment) %>% 
    summarize(
        nCells = n(),
        medTSS = median(TSSEnrichment), 
        min_Frags = min(nFrags),
        mean_Frags = mean(nFrags),
        med_Frags = median(nFrags),
        max_Frags = max(nFrags),
        min_Doub = min(DoubletEnrichment),
        mean_Doub = mean(DoubletEnrichment),
        med_Doub = median(DoubletEnrichment),
        max_Doub = max(DoubletEnrichment),
        quant_Frags_05 = quantile(nFrags, probs = 0.05),
        quant_Frags_10 = quantile(nFrags, probs = 0.10),
        quant_Frags_90 = quantile(nFrags, probs = 0.90),
        quant_Frags_95 = quantile(nFrags, probs = 0.95)) %>% 
    arrange(-med_Frags) %>%
    write.table(
        file = paste0(file.path(outDir.res, qc.table), ".csv"),
        quote = FALSE,
        sep = ",",
        row.names = TRUE,
        col.names = TRUE)

##### Heatmap of gene score to asses quality of cell clusters
markersGS <- getMarkerFeatures(
    ArchRProj = proj_sub, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
proj_sub <- addImputeWeights(proj_sub, reducedDims = "IterativeLSI")

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE)

# ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-QC", width = 8, height = 6, ArchRProj = proj_sub, addDOC = FALSE)

pdf(file.path(outDir.fig, "GeneScores-Marker-Heatmap-QC.pdf"))
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


saveArchRProject(
    ArchRProj = proj_sub, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
