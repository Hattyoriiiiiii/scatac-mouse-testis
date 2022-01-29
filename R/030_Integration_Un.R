suppressPackageStartupMessages({
    library(ArchR)
    library(pheatmap)
    library(tidyverse)
})


#----------------------------------------------------------------------------------------------------
############### Import functions ###############
#----------------------------------------------------------------------------------------------------

source("R/functions.R")
source("R/sceRNA.R")

#----------------------------------------------------------------------------------------------------
############### Set up ###############
#----------------------------------------------------------------------------------------------------

set.seed(42)
addArchRThreads(threads = 32)
addArchRGenome("Mm10")

inDir.ArchRProject <- "processed/02_filtering/01_proj_filtered"
outDir.ArchRProject <- "processed/03_Integration_Un/03_proj_UnInteg"

outDir.proc <- "processed/03_Integration_Un"
outDir.fig <- "figures/03_Integration_Un"
outDir.res <- "results/03_Integration_Un"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

#----------------------------------------------------------------------------------------------------
############### Unconstrained Integration ###############
#----------------------------------------------------------------------------------------------------

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = sceRNA,
    addToArrow = FALSE,
    groupRNA = "AnnotatedClusters",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un")


#----------------------------------------------------------------------------------------------------
############### Assessment of Unconstrained Integration ###############
#----------------------------------------------------------------------------------------------------

cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
cM %>% 
    as.data.frame() %>%
    t() %>%
    write.table(
        file = file.path(outDir.res, "confusionMatrix.txt"), 
        row.names = TRUE, 
        col.names = TRUE, 
        quote = FALSE)

pdf(file.path(outDir.fig, "UMAP-predictedGroup_Un.pdf"))
qc_plot(
    proj,
    qc_metrix = c("Clusters", "predictedGroup_Un", "predictedScore_Un"),
    embedding = "UMAP"
)
dev.off()


saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
