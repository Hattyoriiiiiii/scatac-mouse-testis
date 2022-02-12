suppressPackageStartupMessages({
    library(ArchR)
    library(ComplexHeatmap)
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


inDir.ArchRProject <- "processed/15_trajectory/proj_germ_final"

outDir.proc <- "processed/15_trajectory_subset"
outDir.fig <- "figures/15_trajectory_subset"
outDir.res <- "results/15_trajectory_subset"
# out.table <- "metadata_positiveTFs"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
motifPositions <- getPositions(proj)

cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]
cluster_name <- "ClustersGerm"
embedding <- "UMAP_Germ"

SpgTrajectory <- cellOrder[1:4]
SctTrajectory <- cellOrder[5:9]
StdTrajectory <- cellOrder[10:13]

Spg2SctTrajectory <- cellOrder[3:7]
Sct2StdTrajectory <- cellOrder[7:10]

#----------------------------------------------------------------------------------------------------
############### Add trajectory ###############
#----------------------------------------------------------------------------------------------------

trajName <- c("WholeGermU", "SpgU", "SctU", "StdU", "Spg2SctU", "Sct2StdU")
trajCell <- list(
    cellOrder, SpgTrajectory, SctTrajectory, 
    StdTrajectory, Spg2SctTrajectory, Sct2StdTrajectory)

for (i in 1:length(trajName)) {

    CellTypeU <- trajName[i]

    proj <- addTrajectory(
        ArchRProj = proj, 
        name = trajName[i], 
        groupBy = cluster_name,
        trajectory = trajCell[[i]], 
        embedding = embedding, 
        force = TRUE)

    # positive regulators

    trajMM  <- getTrajectory(
        ArchRProj = proj, 
        name = CellTypeU, 
        useMatrix = "cisBPchromVar_Germ", 
        log2Norm = FALSE)
    
    trajGSM <- getTrajectory(
        ArchRProj = proj, 
        name = CellTypeU, 
        useMatrix = "GeneScoreMatrix", 
        log2Norm = TRUE)

    trajGIM <- getTrajectory(
        ArchRProj = proj, 
        name = CellTypeU, 
        useMatrix = "GeneIntegrationMatrix", 
        log2Norm = FALSE)


    # Positive regulators
    corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
    idx <- grep("^z:", corGIM_MM[[1]]$name2)


    trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1[idx], ]
    trajMM2 <- trajMM[corGIM_MM[[1]]$name2[idx], ]

    trajCombined <- trajGIM2
    combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

    rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))

    ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "coolwarm"),  varCutOff = 0, rowOrder = rowOrder)
    ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)


    pdf(file.path(outDir.fig, paste0("heatmap-", trajName[i], ".pdf")))
    ComplexHeatmap::draw(ht1 + ht2)
    dev.off()


    # embedding, positive tfs, density, gim x mm, gsm x mm
    pdf(file.path(outDir.fig, paste0("UMAP-", trajName[i], ".pdf")))
    p <- plotTrajectory(
        proj, 
        trajectory = trajName[i], 
        embedding = embedding,
        colorBy = "cellColData", 
        continuousSet = "beach",
        name = trajName[i])[[1]]
    print(p)
    dev.off()

    
    ### positive tfs table

    write.table(
        as.vector(corGIM_MM[[1]]$matchname1[idx]), 
        file = file.path(outDir.fig, paste0("tfs-", trajName[i], ".txt")), 
        quote = FALSE, 
        row.names = FALSE, 
        col.names = FALSE)

    ### density ----------

    df <- data.frame(
        celltype = proj@cellColData[[cluster_name]],
        pseudo_time = proj@cellColData[[trajName[i]]]
    ) %>% drop_na() %>% mutate(celltype = factor(celltype, levels = trajCell[[i]]))

    pdf(file.path(outDir.fig, paste0("density-", trajName[i], ".pdf")))
    p <- df %>% 
        ggplot() + 
        geom_density(aes(x=pseudo_time, fill=celltype)) + 
        scale_fill_brewer(palette = "PRGn") +
        theme_classic()
    print(p)
    dev.off()



    pdf(file.path(outDir.fig, paste0("trajectory-", trajName[i], ".pdf")))
    for (motifs in as.vector(corGIM_MM[[1]]$matchname1[idx])) {

        markerMotifs <- unlist(lapply(paste0("^", motifs, "_"), function(x) grep(x, names(motifPositions), value = TRUE)))

        p1 <- plotTrajectory(
            proj, 
            trajectory = trajName[i], 
            embedding = embedding,
            colorBy = "GeneIntegrationMatrix", 
            name = motifs,
            imputeWeights = getImputeWeights(proj),
            addArrow = FALSE,
            continuousSet = "coolwarm")
        p2 <- plotTrajectory(
            proj, 
            trajectory = trajName[i], 
            embedding = embedding,
            colorBy = "cisBPchromVar_Germ",
            name = paste0("z:", markerMotifs),
            imputeWeights = getImputeWeights(proj),
            addArrow = FALSE,
            continuousSet = "solarExtra")

        p <- ggAlignPlots(p1[[2]], p2[[2]], type="h")
        print(p)
    }

    dev.off()
}
