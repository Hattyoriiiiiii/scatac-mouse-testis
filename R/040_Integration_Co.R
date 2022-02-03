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

inDir.ArchRProject <- "processed/03_Integration_Un/03_proj_UnInteg"
outDir.ArchRProject <- "processed/04_Integration_Co/04_proj_CoInteg"

outDir.proc <- "processed/04_Integration_Co"
outDir.fig <- "figures/04_Integration_Co"
outDir.res <- "results/04_Integration_Co"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

#----------------------------------------------------------------------------------------------------
############### Constrained Integration ###############
#----------------------------------------------------------------------------------------------------

cSertoli <- c("Sertoli", "Sertoli_1", "Sertoli_2")
cSomatic <- c("PTM", "Endothelial_cells", "Leydig_1", "Leydig_2", "Fetal_Leydig_1", "Fetal_Leydig_2", "Immature_Leydig", "Interstitial_tMg")
cSpg2Sct <- c("Undifferentiated_Spermatogonia", "Spermatogonia", "Differentiating_Spermatogonia_1", "Differentiating_Spermatogonia_2", "Leptotene", "Zygotene", "eP1", "eP2")
cLateScyte <- c("mP", "lP1", "lP2", "D", "MI", "MII")
cStid <- c(paste0("S", 1:11))

clustSertoli <- c("C8")
clustSomatic <- c("C9")
clustSpg2Sct <- c("C3", "C4", "C5", "C6", "C7")
clustLateScyte <- c("C2")
clustStid <- c("C10", "C11", "C12", "C13")

rnaSertoli <- colnames(sceRNA)[sceRNA$AnnotatedClusters %in% cSertoli]
rnaSomatic <- colnames(sceRNA)[sceRNA$AnnotatedClusters %in% cSomatic]
rnaSpg2Sct <- colnames(sceRNA)[sceRNA$AnnotatedClusters %in% cSpg2Sct]
rnaLateScyte <- colnames(sceRNA)[sceRNA$AnnotatedClusters %in% cLateScyte]
rnaStid <- colnames(sceRNA)[sceRNA$AnnotatedClusters %in% cStid]

groupList <- SimpleList(
    SC = SimpleList(
        ATAC = proj$cellNames[proj$Clusters %in% clustSertoli],
        RNA = rnaSertoli
    ),
    Somatic = SimpleList(
        ATAC = proj$cellNames[proj$Clusters %in% clustSomatic],
        RNA = rnaSomatic
    ),
    Spg2Sct = SimpleList(
        ATAC = proj$cellNames[proj$Clusters %in% clustSpg2Sct],
        RNA = rnaSpg2Sct
    ),
    lScyte = SimpleList(
        ATAC = proj$cellNames[proj$Clusters %in% clustLateScyte],
        RNA = rnaLateScyte
    ),
    Stid = SimpleList(
        ATAC = proj$cellNames[proj$Clusters %in% clustStid],
        RNA = rnaStid
    )
)


proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = sceRNA,
    addToArrow = FALSE,
    force = TRUE,
    groupList = groupList,
    groupRNA = "AnnotatedClusters",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co")


#----------------------------------------------------------------------------------------------------
############### Assessment of Constrained Integration ###############
#----------------------------------------------------------------------------------------------------

cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Co))
cM %>% 
    as.data.frame() %>%
    t() %>%
    write.table(
        file = file.path(outDir.res, "confusionMatrix.txt"), 
        row.names = TRUE, 
        col.names = TRUE, 
        quote = FALSE)

pdf(file.path(outDir.fig, "UMAP-predictedGroup_Co.pdf"))
qc_plot(
    proj,
    qc_metrix = c("Clusters", "predictedGroup_Co", "predictedScore_Co"),
    embedding = "UMAP"
)
dev.off()


proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = sceRNA,
    addToArrow = TRUE,
    force = TRUE,
    groupList = groupList,
    groupRNA = "AnnotatedClusters",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore")

proj <- addImputeWeights(proj)

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
