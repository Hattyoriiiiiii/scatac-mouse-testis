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

inDir.ArchRProject <- "processed/04_Integration_Co/04_proj_CoInteg"
outDir.ArchRProject <- "processed/05_peakcall/05_proj_peakcall"

outDir.proc <- "processed/05_peakcall"
outDir.fig <- "figures/plots"

mkdir(outDir.proc)
mkdir(outDir.fig)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

markerGenes <- read.table("gene_sets/RepresentativeGenes.txt")[["V1"]]
cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cluster_name <- "ClustersWhole"

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

remapClust <- c(
    "C2" = "mlP_D", 
    "C3" = "Z_eP", 
    "C4" = "preLep",
    "C5" = "Lep", 
    "C6" = "Diff",
    "C7" = "Undiff",
    "C8" = "Sertoli",
    "C9" = "Soma",
    "C10" = "ES", 
    "C11" = "lRS", 
    "C12" = "eRS", 
    "C13" = "mRS"
)

labelNew <- mapLabels(names(remapClust), oldLabels = names(remapClust), newLabels = remapClust)
proj$ClustersWhole <- mapLabels(
    proj$Clusters, 
    newLabels = labelNew, 
    oldLabels = names(remapClust))

#----------------------------------------------------------------------------------------------------
############### Dotplot of marker genes ###############
#----------------------------------------------------------------------------------------------------

gs_mat <- getMatrixHeatmap(
    ArchRProj = proj,
    ClustersName = cluster_name,
    MatrixName = "GeneScoreMatrix",
    markerGenes = markerGenes,
    cellOrder = cellOrder,
    motif = FALSE
)


pdf(file.path(outDir.fig, "dotplot_markergenes.pdf"))
gs_mat[markerGenes, rev(cellOrder)] %>% 
    t() %>% 
    as.matrix() %>% 
    reshape2::melt() %>% 
    ggplot(aes(x = Var2, y = Var1, color = value, size = value)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=1)) +
    theme(axis.ticks = element_blank()) +
    xlab('') + ylab('') + scale_size(range = c(0,10)) +
    scale_color_continuous(type = "viridis") +
    guides(
        color= guide_legend(), 
        size=guide_legend())
dev.off()

#----------------------------------------------------------------------------------------------------
############### Calling peaks ###############
#----------------------------------------------------------------------------------------------------

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersWhole")

pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "ClustersWhole", 
    pathToMacs2 = pathToMacs2)

proj <- addPeakMatrix(proj)
saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)