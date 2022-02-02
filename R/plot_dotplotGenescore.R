suppressPackageStartupMessages({
    library(ArchR)
    library(pheatmap)
    library(tidyverse)
    library(ComplexHeatmap)
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

outDir.fig <- "figures/plots"
mkdir(outDir.fig)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

markerGenes <- read.table("gene_sets/positiveTFs_whole.txt")[["V1"]]
markerGenes <- c(
    "Nanos2",
    "Sall4",
    "Stra8",
    "Prdm9",
    "Spo11",
    "Dmc1",
    "Prss50",
    "Piwil1",
    "Adam3",
    "Pgk2",
    "Acrv1",
    "Prm2",
    "Sox9",
    "Cldn11",
    "Fabp3",
    "Inhba",
    "Sycp3",
    "Rhox13",
    "Meiob",
    "Gata4",
    "Neurog3",
    "Lin28a"
)

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cluster_name <- "ClustersWhole"

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