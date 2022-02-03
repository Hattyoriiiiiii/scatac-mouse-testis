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

outDir.proc <- "processed/01_qc"
outDir.fig <- "figures/01_qc"
outDir.res <- "results/01_qc"
outDir.ArchRProject <- "processed/01_qc/01_proj_qc"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

qc.table <- "QC_metadata_before_filtering"

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

sample.name <- c("Mouse_4weeks", "Mouse_7weeks")
inputFiles <- c("processed/scATAC/Mouse_4weeks.tsv.gz", 
                "processed/scATAC/Mouse_7weeks.tsv.gz")

names(inputFiles) <- sample.name

markerGenes <- read.table("gene_sets/markerGenes.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### 1 st filtering ###############
#----------------------------------------------------------------------------------------------------

# Creating An ArrowFiles
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    filterTSS = 4,  # minTSS
    filterFrags = 1000,  # minFrags
    QCDir = outDir.res,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE)

# Inferring doublets
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "UMAP",
    LSIMethod = 1,
    outDir = outDir.res)

# Creating An ArchRProject
proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = outDir.ArchRProject,
    copyArrows = TRUE, 
    showLogo = FALSE)

# removing top 5% of cells using filterRatio = 1
proj <- filterDoublets(proj, filterRatio = 1)

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1, p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, 
        width = 4, height = 4)

pdf(file.path(outDir.fig, "QC-Sample-FragSizes-TSSProfile.pdf"))
print(ggAlignPlots(p1, p2, type = "v"))
dev.off()

#----------------------------------------------------------------------------------------------------
############### 2 nd filtering ###############
#----------------------------------------------------------------------------------------------------

############### resolution : 2.0 ###############
### Purpose
### 目的はdoubletsやlow qualiy cellsがエンリッチしているクラスターを除去することである。
### そのため、比較的高いresolutionでSeurat's Louvain clusteringを行う

# Iterative Latent Semantic Indexing
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI_qc", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1, 0.2, 0.4, 0.8), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    excludeChr = c("chrX", "chrY", "chrMT"),
    force = TRUE)

# Clustering using Seurat’s FindClusters() function
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI_qc",
    method = "Seurat",
    name = "Clusters_qc",
    resolution = 2.0,
    force = TRUE)

# UMAPアルゴリズムを用いて可視化することでputative doublets / low quality cells の分布を把握する
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI_qc", 
    name = "UMAP_qc", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE)


#----------------------------------------------------------------------------------------------------
############### Quality Assessment ###############
#----------------------------------------------------------------------------------------------------
# 1. Sample : batch effectが存在するかの確認。
# 2. Cluster : クラスターと番号の関係を確認。これは上のdataframeと合わせて確認。
# 3. Fragment数 : fragment数が多い細胞はdoubletである可能性が高い。ただし、spermatidsなど少なくなる可能性がある。
# 4. DoubletEnrichment : Doubletである可能性が高い。
# 5. TSS enrichment : low qualiy である可能性。ただしGerm cellsの場合、そうとは限らない。
#----------------------------------------------------------------------------------------------------

pdf(file.path(outDir.fig, "QC-second-filtering.pdf"))
qc_plot(proj, embedding = "UMAP_qc")
dev.off()


pdf(file.path(outDir.fig, "QC-before-subsetting.pdf"))
p <- ggPoint(
        x=proj$nFrags,
        y=proj$DoubletEnrichment,
        colorDensity = T,
        rastr = T) +
    geom_hline(yintercept = 3, col="red", lty="dashed") +
    geom_vline(xintercept = 30000, col="red", lty="dashed") +
    xlab("nFragments") +
    ylab("DoubletEnrichment")
print(p)
dev.off()

# save the quality metrix
getCellColData(proj) %>%
    as.data.frame() %>% 
    group_by(Clusters_qc) %>% 
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
        row.names = FALSE,
        col.names = TRUE)

##### Heatmap of gene score to asses quality of cell clusters
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters_qc",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
proj <- addImputeWeights(proj, reducedDims = "IterativeLSI_qc")

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE)

# ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-QC", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

pdf(file.path(outDir.fig, "GeneScores-Marker-Heatmap-QC.pdf"))
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)