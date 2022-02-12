suppressPackageStartupMessages({
    library(ArchR)
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

inDir.ArchRProject <- "processed/20_subsetGerm_peakAnno/20_proj_germ_peakAnno"
outDir.ArchRProject <- "processed/30_subsetGerm_p2g/07_proj_germ_p2g"

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)
cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]

outDir.proc <- "processed/30_subsetGerm_p2g"
mkdir(outDir.proc)

#----------------------------------------------------------------------------------------------------
############### get peak2gene for germ cells ###############
#----------------------------------------------------------------------------------------------------

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI_Germ")

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
