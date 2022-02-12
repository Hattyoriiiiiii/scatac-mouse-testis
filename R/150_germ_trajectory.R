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

inDir.ArchRProject <- "processed/30_subsetGerm_p2g/07_proj_germ_p2g"
outDir.ArchRProject <- "processed/15_trajectory/proj_germ_final"

outDir.proc <- "processed/15_trajectory"
outDir.fig <- "figures/15_trajectory"
outDir.res <- "results/15_trajectory"
# out.table <- "metadata_positiveTFs"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]
GermTrajectory <- cellOrder[1:13]
cluster_name <- "ClustersGerm"

#----------------------------------------------------------------------------------------------------
############### Add trajectory ###############
#----------------------------------------------------------------------------------------------------

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "GermU", 
    groupBy = cluster_name,
    trajectory = GermTrajectory, 
    embedding = "UMAP_Germ", 
    force = TRUE)

p <- plotTrajectory(
    proj, 
    trajectory = "GermU", 
    colorBy = "cellColData", 
    name = "GermU"
)

plotPDF(
    p, 
    name = "Plot-Germ-Traj-UMAP.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, 
    width = 5, height = 5)


#----------------------------------------------------------------------------------------------------
############### density along with pseudotime ###############
#----------------------------------------------------------------------------------------------------

celltype <- proj@cellColData[[cluster_name]]
time <- proj@cellColData$GermU

df <- data.frame(
    celltype = celltype,
    pseudo_time = time
) %>% mutate(celltype = factor(celltype, levels = GermTrajectory))

p_celltype <- df %>% 
  ggplot() + 
  geom_density(aes(x=pseudo_time, fill=celltype)) + 
  scale_fill_brewer(palette = as.vector(ArchRPalettes$stallion2)[1:13]) +
  theme_classic()

pdf(file.path(outDir.fig, "density_pseudotime.pdf"))
print(p_celltype); print(p)
dev.off()

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)