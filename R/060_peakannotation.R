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

inDir.ArchRProject <- "processed/05_peakcall/05_proj_peakcall"
outDir.ArchRProject <- "processed/06_peakannotation/06_proj_peakannotation"

outDir.proc <- "processed/06_peakannotation"
outDir.fig <- "figures/06_peakannotation"
outDir.res <- "results/06_peakannotation"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "ClustersWhole",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

saveRDS(markersPeaks, file.path(outDir.proc, "se_markerPeaks.rds"))

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks, 
    nPrint = 0,
    nLabel = 0,
    cutOff = "FDR <= 0.1 & Log2FC >= 1",
    transpose = TRUE)

pdf(file.path(outDir.fig, "Peak-Marker-Heatmap.pdf"))
draw(heatmapPeaks, 
    heatmap_legend_side = "bot", 
    annotation_legend_side = "bot",
    row_order = cellOrder)
dev.off()

plotPDF(heatmapPeaks, 
        name = "Peak-Marker-Heatmap", 
        width = 8, height = 6, 
        ArchRProj = proj, 
        addDOC = FALSE)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- addMotifAnnotations(
    ArchRProj = proj, 
    motifSet = "cisbp",
    name = "Motif_cisbp", 
    force = TRUE)

proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Motif_cisbp",
    matrixName = "cisBPchromVar",
    force = TRUE)

# motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif_cisbp",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

pdf(file.path(outDir.fig, "Motifs-Enriched-Marker-Heatmap.pdf"))
ComplexHeatmap::draw(heatmapEM, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot", 
                     row_order = cellOrder)
dev.off()

plotPDF(heatmapEM, 
        name = "Motifs-Enriched-Marker-Heatmap", 
        width = 8, height = 6, 
        ArchRProj = proj, 
        addDOC = FALSE)

# #----------------------------------------------------------------------------------------------------
# ############### TFBS ###############
# #----------------------------------------------------------------------------------------------------

# proj <- addArchRAnnotations(ArchRProj = proj, collection = "EncodeTFBS")
# enrichEncode <- peakAnnoEnrichment(
#     seMarker = markersPeaks,
#     ArchRProj = proj,
#     peakAnnotation = "EncodeTFBS",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#   )

# heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)

# pdf(file.path(outDir.fig, "EncodeTFBS-Enriched-Marker-Heatmap.pdf"))
# ComplexHeatmap::draw(heatmapEncode, 
#     heatmap_legend_side = "bot", 
#     annotation_legend_side = "bot",
#     row_order = cellOrder)
# dev.off()

# plotPDF(heatmapEncode, 
#         name = "EncodeTFBS-Enriched-Marker-Heatmap", 
#         width = 8, height = 6, 
#         ArchRProj = proj, 
#         addDOC = FALSE)


saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)