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

outDir.proc <- "processed/14_positiveRegulators"
outDir.fig <- "figures/14_positiveRegulators"
outDir.res <- "results/14_positiveRegulators"
out.table <- "metadata_positiveTFs"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_germ.txt")[["V1"]]
cluster_name <- "ClustersGerm"

#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

seGroupMotif <- getGroupSE(
    ArchRProj = proj, 
    useMatrix = "cisBPchromVar_Germ", 
    groupBy = cluster_name)

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
        rowMaxs(assay(seZ) - assay(seZ)[,x])
    }) %>% Reduce("cbind", .) %>% rowMaxs

corGIM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "cisBPchromVar_Germ",
    reducedDims = "IterativeLSI_Germ")

# GeneIntegration
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$cisBPchromVar_Germ_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"cisBPchromVar_Germ_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75, , na.rm = TRUE))] <- "YES"

corGIM_MM %>%
    write.table(
        file = paste0(file.path(outDir.res, out.table), ".csv"),
        quote = FALSE,
        sep = ",",
        row.names = FALSE,
        col.names = TRUE)

df <- read.csv(paste0(file.path(outDir.res, out.table), ".csv"))
write.table(
    df[df$TFRegulator == "YES",]$GeneIntegrationMatrix_name, 
    file="gene_sets/positiveTFs_germ.txt", 
    quote=FALSE, 
    row.names=FALSE, 
    col.names=FALSE)

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

pdf(file.path(outDir.fig, "PositiveTFs_Germ.pdf"))
print(p)
dev.off()