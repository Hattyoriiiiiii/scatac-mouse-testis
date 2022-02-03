suppressPackageStartupMessages({
    library(ArchR)
    library(ChIPpeakAnno)
    library(ChIPseeker)
    library(clusterProfiler)
    library(ComplexHeatmap)
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(org.Mm.eg.db)
    library(pheatmap)
    library(RColorBrewer)
    library(ReactomePA)
    library(rtracklayer)
    library(stringr)
    library(tidyr)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
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

outDir.proc <- "processed/07_p2g_distal"
outDir.fig <- "figures/07_p2g_distal"
outDir.res <- "results/07_p2g_distal"

mkdir(outDir.proc)
mkdir(outDir.fig)
mkdir(outDir.res)


##### parameters ---------------
km_number <- 11
EXP <- "exp_01"
##### --------------------------

PARDIR <- file.path(outDir.res, sprintf("CREs/%s", EXP))
DIRList <- list("bed", "go", "metaplot", "motif")

if (!dir.exists(PARDIR)) {
    for (DIR in DIRList) {
        dir.create(file.path(PARDIR, DIR), recursive=TRUE)
        # print(file.path(PARDIR, DIR))
    }
}

#----------------------------------------------------------------------------------------------------
############### Inputs ###############
#----------------------------------------------------------------------------------------------------

proj <- loadArchRProject(path = inDir.ArchRProject, showLogo = FALSE)

cellOrder <- read.table("gene_sets/cellOrder_whole.txt")[["V1"]]

color <- as.vector(ArchRPalettes$stallion2)[1:14]
names(color) <- cellOrder
my_color <- list(Cell.Type = color)


#----------------------------------------------------------------------------------------------------
############### Cell-type definition based on gene score and signal track ###############
#----------------------------------------------------------------------------------------------------

# Peak2GeneLinksを取り出す
p2g <- plotPeak2GeneHeatmap(
    ArchRProj = proj, 
    corCutOff = 0.45,
    groupBy = "ClustersWhole",
    k = km_number,
    returnMatrices = TRUE)


p2g_RNA <- p2g$RNA$matrix
p2g_ATAC <- p2g$ATAC$matrix


# ここでp2g_idのfilteringを行う
df <- p2g$Peak2GeneLinks %>% 
    as.data.frame() %>% 
    dplyr::select(peak, gene, Correlation)

df$p2g_id <- rownames(p2g$Peak2GeneLinks)  # idの付加

df[c("chr", "start_end")] <- str_split_fixed(df$peak, ":", 2)  # chr, start, endの分割
df[c("start", "end")] <- str_split_fixed(df$start_end, "-", 2) # chr, start, endの分割
df <- GRanges(df[!duplicated(df$peak),])

peakSet <- unique(proj@peakSet)
peakSet$CellType <- names(peakSet)
names(peakSet) <- NULL
peakSet$uniq_id <- paste0(peakSet@seqnames, ":", peakSet@ranges)
peakSet <- peakSet[!duplicated(peakSet$uniq_id, fromLast = TRUE),]
ol <- findOverlapsOfPeaks(df, peakSet, ignore.strand=FALSE)
metadata_p2g <- ol$overlappingPeaks[[1]]

p2g_distal <- metadata_p2g[metadata_p2g$peakType %in% c("Distal", "Intronic", "Exonic"),]
p2g_proximal <- metadata_p2g[metadata_p2g$peakType %in% c("Promoter"),]



BroadAnno <- list(
    SPG = cellOrder[1:2],
    SCT = cellOrder[3:6],
    STD = cellOrder[7:10],
    SOMA = cellOrder[11:12])

remap <- function(cell) {
    if (cell %in% BroadAnno$SPG) {
        celltype <- "SPG"
    } else if (cell %in% BroadAnno$SCT) {
        celltype <-  "SCT"
    } else if (cell %in% BroadAnno$STD) {
        celltype <- "STD"  
    } else if (cell %in% BroadAnno$SOMA) {
        celltype <- "SOMA"  
    }
    return(celltype)
}

# 初期状態のmatrixのcolumn idx
annot_cols <- data.frame(
    "Cell Type" = p2g$RNA$colData$groupBy,
    row.names = rownames(p2g$RNA$colData))
annot_cols$Cell.Type <- factor(annot_cols$Cell.Type, levels = cellOrder)
annot_cols <- annot_cols %>% arrange(annot_cols$Cell.Type)
annot_cols$BroadAnno <- lapply(annot_cols$Cell.Type, FUN = remap)
annot_cols$BroadAnno <- factor(annot_cols$BroadAnno, levels = names(BroadAnno))


min_val <- -2
max_val <- 2


##### ---------------
# 後にソートするため、クラスターごとのidを取り出してリストとして保持する

annot_rows_list <- list()
RNA_df_list <- list()
ATAC_df_list <- list()
cres_gr_list <- list()

for (i in 1:km_number) {
    p2g_id <- p2g$Peak2GeneLinks[p2g$RNA$kmeansId == i, ] %>% rownames()
    
    # ここでfilteringする
    # p2g_dfでfilteringした後のp2g_idを用いる
    ### correlationの正負 -> Enhancer/Silencer
    ### Distal/Proximal
    # p2g_id[p2g_id %in% p2g_df$p2g_id]
    p2g_id <- p2g_id[p2g_id %in% p2g_distal$p2g_id]
    
    # ----- Export of bed files
    
    file.name <- file.path(outDir.res, sprintf("CREs/%s/bed/cluster_%02d.bed", EXP, i))
    gr_p2g <- p2g$Peak2GeneLinks[p2g_id, ]

    chr <- lapply(gr_p2g$peak, function(x) {strsplit(x, ":")[[1]][1]}) %>% as.character() %>% as.vector()
    range <- lapply(gr_p2g$peak, function(x) {strsplit(x, ":")[[1]][2]})
    start <- lapply(range, function(x) {strsplit(x, "-")[[1]][1]}) %>% as.numeric() %>% as.vector()
    end <- lapply(range, function(x) {strsplit(x, "-")[[1]][2]}) %>% as.numeric() %>% as.vector()

    gr <- GRanges(seqnames = chr, ranges = IRanges(start, end))
    gr$kmeans <- i
    gr$gene <- gr_p2g$gene
    gr$FDR <- gr_p2g$FDR
    gr$Correlation <- gr_p2g$Correlation
#     gr_list[[sprintf("gr_%d", i)]] <- unique(gr)
    
    # duplicateを除きたい
    export.bed(sort(gr), file.name)
    
#     cres.name <- sprintf("cluster_%02d", i)
#     cres_gr_list[[cres.name]] <- gr
    cres_gr_list[[i]] <- gr
    
    # ----- 
    
    file.go.name.svg <- file.path(outDir.res, paste0("CREs/", EXP, sprintf("/go/go_cluster_%02d.svg", i)))
    file.go.name.png <- file.path(outDir.res, paste0("CREs/", EXP, sprintf("/go/go_cluster_%02d.png", i)))
    file.cnet.name.svg <- file.path(outDir.res, paste0("CREs/", EXP, sprintf("/go/cnet_cluster_%02d.svg", i)))
    file.cnet.name.png <- file.path(outDir.res, paste0("CREs/", EXP, sprintf("/go/cnet_cluster_%02d.png", i)))
    print(i)

    symbols <- p2g$Peak2GeneLinks[p2g$ATAC$kmeansId == i, ]$gene
    FDR <- p2g$Peak2GeneLinks[p2g$ATAC$kmeansId == i, ]$FDR

    genes <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = as.character(symbols),
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "SYMBOL")

    genes.ids <- genes[, "ENTREZID"]


    # こっから先は関数化 ------------------------------
    
    ########## enrichedGO ##########
    enriched_go <- clusterProfiler::enrichGO(
        gene = genes.ids,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        readable = TRUE)

    enriched_go.simple <- simplify(enriched_go)
    # return (enriched_go.simple)


    ########## plot関数 ##########
    
    ### dotplot ----------
    p_go <- clusterProfiler::dotplot(
        enriched_go.simple
    ) + scale_y_discrete(labels=function(x) str_wrap(x, width = 25))
    ggsave(plot = p_go, filename = file.go.name.svg)
    ggsave(plot = p_go, filename = file.go.name.png)
#     print(p_go)

    ### cnetplot ----------
    p_cnet <- clusterProfiler::cnetplot(
        enriched_go.simple, 
        categorySize="pvalue",
        node_label = "all",
        cex_label_category = 1.4)
    ggsave(plot = p_cnet, filename = file.cnet.name.svg)
    ggsave(plot = p_cnet, filename = file.cnet.name.png)
#     print(p_cnet)
    
#      # こっから前は関数化 ------------------------------

    # クラスターごとにmatrixを取り出す
    RNA_mat <- p2g_RNA[p2g_id, ]
    ATAC_mat <- p2g_ATAC[p2g_id, ]
    
    # capping
    RNA_mat[which(RNA_mat < min_val)] = min_val
    RNA_mat[which(RNA_mat > max_val)] = max_val
    ATAC_mat[which(ATAC_mat < min_val)] = min_val
    ATAC_mat[which(ATAC_mat > max_val)] = max_val
    
    annot_rows <- data.frame(
            row.names = p2g_id)
    annot_rows$cluster <- as.character(i)
#     annot_rows$cluster <- as.character(i)
    
    ### --------------- 
    annot_rows_list[[i]] <- annot_rows
    
    # columnの順番を分化進行
    RNA_df_list[[i]] <- RNA_mat[, rownames(annot_cols)] %>% as.data.frame()
    ATAC_df_list[[i]] <- ATAC_mat[, rownames(annot_cols)] %>% as.data.frame()
}
##### ---------------




##### ---------------
# heatmapのクラスターを、最大値が分化過程の初期であるほど上に来るようにソートする
clustMax <- c()  # idx
for (i in 1:km_number) {
    c <- colMeans(ATAC_df_list[c(i)][[1]])
    c <- which.max(c) %>% as.numeric()
    clustMax <- c(clustMax, c)
}

# matrixの行の順番をソートする
RNA <- bind_rows(RNA_df_list[sort.list(clustMax)]) %>% as.matrix()
ATAC <- bind_rows(ATAC_df_list[sort.list(clustMax)]) %>% as.matrix()

cres.grl <- GRangesList(cres_gr_list[sort.list(clustMax)])
saveRDS(cres.grl, file.path(outDir.res, "cres.grl.rds"))

##### ---------------

annot_rows <- bind_rows(annot_rows_list[sort.list(clustMax)])
row_level <- annot_rows$cluster %>% unique() %>% as.character()
annot_rows$cluster <- factor(annot_rows$cluster, level = row_level)



##### Plot heatmap ---------

pdf(file.path(outDir.res, sprintf("CREs/%s/heatmap.pdf", EXP)))
split_rows = factor(annot_rows$cluster)
split_cols = factor(annot_cols$BroadAnno)

ha_rna = HeatmapAnnotation(BroadAnnotation = annot_cols$BroadAnno,
                       CellType = annot_cols$Cell.Type, 
                       border = c(CellType = TRUE, BroadAnnotation = TRUE),
                       show_annotation_name = FALSE)

ha_atac = HeatmapAnnotation(BroadAnnotation = annot_cols$BroadAnno,
                       CellType = annot_cols$Cell.Type, 
                       border = c(CellType = TRUE, BroadAnnotation = TRUE),
                       show_annotation_name = TRUE)

ht_rna <- Heatmap(RNA, name = "RNA", 
        row_split = split_rows, 
        column_split = split_cols,
        cluster_row_slices = FALSE,
        border = "#404040",
        top_annotation = ha_rna, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        column_title = "RNA")

ht_atac <- Heatmap(ATAC, name = "ATAC", 
        row_split = split_rows, 
        column_split = split_cols,
        border = "#404040",
        top_annotation = ha_atac, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        column_title = "ATAC")

ht_rna + ht_atac

dev.off()