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

inDir.ArchRProject <- "processed/06_peakannotation/06_proj_peakannotation"
outDir.ArchRProject <- "processed/07_p2g/07_proj_p2g"

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI")

saveArchRProject(
    ArchRProj = proj, 
    outputDirectory = outDir.ArchRProject, 
    load = FALSE)
