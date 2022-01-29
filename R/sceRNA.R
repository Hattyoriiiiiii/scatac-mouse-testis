library(SingleCellExperiment)

#----------------------------------------------------------------------------------------------------
############### Create SingleCellExperiment obj ###############
#----------------------------------------------------------------------------------------------------

sceRNA <- readRDS("processed/scRNA/SCE_emptyDrops.rds")
rownames(sceRNA) <- rowData(sceRNA)$Symbol

filter <- (!is.na(sceRNA@colData$P15Clusters)) & sceRNA@colData$AnnotatedClusters == "NewCell"
sceRNA@colData[filter, ]$AnnotatedClusters <- sceRNA@colData[filter, ]$P15Clusters