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

inDir.ArchRProject <- "processed/07_p2g/07_proj_p2g"

outDir.proc <- "processed/08_positiveRegulators"
outDir.fig <- "figures/08_positiveRegulators"
outDir.res <- "results/08_positiveRegulators"
out.table <- "metadata_positiveTFs"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

seGroupMotif <- getGroupSE(
    ArchRProj = proj, 
    useMatrix = "cisBPchromVar", 
    groupBy = "ClustersWhole")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
        rowMaxs(assay(seZ) - assay(seZ)[,x])
    }) %>% Reduce("cbind", .) %>% rowMaxs

corGIM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "cisBPchromVar",
    reducedDims = "IterativeLSI")

# GeneIntegration
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$cisBPchromVar_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"cisBPchromVar_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75, , na.rm = TRUE))] <- "YES"

corGIM_MM %>%
    write.table(
        file = paste0(file.path(outDir.res, out.table), ".csv"),
        quote = FALSE,
        sep = ",",
        row.names = TRUE,
        col.names = TRUE)


p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.1)
  ) + 
  ggrepel::geom_label_repel(
      data = data.frame(corGIM_MM)[data.frame(corGIM_MM)$TFRegulator == "YES",], 
      aes(x=cor, y=maxDelta, label=GeneIntegrationMatrix_name),
      size = 2.5,
      nudge_x = 1,
      color = "black",
      max.overlaps = 100
  )

pdf(file.path(outDir.fig, "PositiveTFs.pdf"))
print(p)
dev.off()