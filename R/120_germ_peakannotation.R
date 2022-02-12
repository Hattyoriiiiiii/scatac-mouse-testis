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

inDir.ArchRProject <- "processed/20_subsetGerm_peak/20_proj_germ_peak"
outDir.ArchRProject <- "processed/20_subsetGerm_peakAnno/20_proj_germ_peakAnno"

outDir.proc <- "processed/20_subsetGerm_peakAnno"
outDir.fig <- "figures/20_subsetGerm_peakAnno"
outDir.res <- "results/20_subsetGerm_peakAnno"
out.table <- "metadata_germ"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Plotting ###############
#----------------------------------------------------------------------------------------------------

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "ClustersGerm",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

saveRDS(markersPeaks, file.path(outDir.proc, "se_markerPeaks.rds"))

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks, 
    nPrint = 0,
    nLabel = 0,
    cutOff = "FDR <= 0.1 & Log2FC >= 1",
    transpose = TRUE)

pdf(file.path(outDir.fig, "Peak-Marker-Heatmap-Germ.pdf"))
draw(heatmapPeaks, 
    heatmap_legend_side = "bot", 
    annotation_legend_side = "bot",
    row_order = cellOrder)
dev.off()

plotPDF(heatmapPeaks, 
        name = "Peak-Marker-Heatmap-Germ", 
        width = 8, height = 6, 
        ArchRProj = proj, 
        addDOC = FALSE)


#----------------------------------------------------------------------------------------------------
############### Peak Annotation ###############
#----------------------------------------------------------------------------------------------------

proj <- addMotifAnnotations(
    ArchRProj = proj, 
    motifSet = "cisbp",
    name = "Motif_cisbp_Germ", 
    force = TRUE)

proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Motif_cisbp_Germ",
    matrixName = "cisBPchromVar_Germ",
    force = TRUE)

# motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif_cisbp_Germ",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 3, transpose = TRUE)

pdf(file.path(outDir.fig, "Motifs-Enriched-Marker-Heatmap-Germ.pdf"))
ComplexHeatmap::draw(heatmapEM, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot", 
                     row_order = cellOrder)
dev.off()

plotPDF(heatmapEM, 
        name = "Motifs-Enriched-Marker-Heatmap-Germ", 
        width = 8, height = 6, 
        ArchRProj = proj, 
        addDOC = FALSE)

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
