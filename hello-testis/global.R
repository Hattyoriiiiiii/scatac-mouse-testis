suppressPackageStartupMessages({
    library(ArchR)
})

#----------------------------------------------------------------------------------------------------
############### Set up ###############
#----------------------------------------------------------------------------------------------------

set.seed(42)
addArchRThreads(threads = 32)
addArchRGenome("Mm10")

inDir.ArchRProject <- "/work/proj_germ_final"

print("Please wait...")
proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)


motifPositions <- getPositions(proj)
genes <- readRDS("/work/hello-testis/genes.rds")
cellOrder <- read.table("/work/gene_sets/cellOrder_germ.txt")[["V1"]]
cluster_name <- "ClustersGerm"
embedding <- "UMAP_Germ"
trajectory <- "WholeGermU"

SpgTrajectory <- cellOrder[1:4]
SctTrajectory <- cellOrder[5:9]
StdTrajectory <- cellOrder[10:13]
Spg2SctTrajectory <- cellOrder[3:7]
Sct2StdTrajectory <- cellOrder[7:10]

i <- 1
trajName <- c("WholeGermU", "SpgU", "SctU", "StdU", "Spg2SctU", "Sct2StdU")
trajCell <- list(
    cellOrder, SpgTrajectory, SctTrajectory, 
    StdTrajectory, Spg2SctTrajectory, Sct2StdTrajectory)

proj <- addTrajectory(
    ArchRProj = proj, 
    name = trajName[i], 
    groupBy = cluster_name,
    trajectory = trajCell[[i]], 
    embedding = embedding, 
    force = TRUE)

peakSet <- proj@peakSet
TFBSs <- read.table("/work/hello-testis/TFBSlists.txt")[["V1"]]
