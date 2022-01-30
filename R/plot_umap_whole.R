suppressPackageStartupMessages({
    library(ArchR)
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

inDir.ArchRProject <- "processed/06_peakannotation/06_proj_peakannotation"

outDir.fig <- "figures/plots"
mkdir(outDir.fig)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]
cluster_name <- "ClustersWhole"

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------


Barcode_UMAP <- proj@embeddings$UMAP$df %>% 
    tibble::rownames_to_column("Barcodes")
colnames(Barcode_UMAP) <- c("Barcodes", "UMAP_1", "UMAP_2")

metadata <- proj@cellColData %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Barcodes")

cluster_length <- metadata[[cluster_name]] %>%
    unique() %>% 
    length()


pdf(file.path(outDir.fig, "UMAP_whole.pdf"))

p <- ggplot(metadata %>% left_join(Barcode_UMAP, by = "Barcodes"), aes(x = UMAP_1, y = UMAP_2, color = ClustersWhole)) + 
    geom_point(size = 0.5) + 
    coord_fixed(ratio = 1.4) +
    theme_classic() + 
    theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
    guides(colour=guide_legend(override.aes = list(alpha=1, size=7, shape=15))) +
    scale_color_manual(values=as.vector(ArchRPalettes$stallion2)[1:cluster_length], limits = cellOrder)
print(p)
dev.off()