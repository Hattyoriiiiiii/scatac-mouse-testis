
#----------------------------------------------------------------------------------------------------
############### Import functions ###############
#----------------------------------------------------------------------------------------------------

mkdir <- function(path) {
    if (!dir.exists(path)) {
       dir.create(path, recursive = TRUE)
    }
}

#----------------------------------------------------------------------------------------------------
############### Plotters ###############
#----------------------------------------------------------------------------------------------------

qc_plot <- function(
    ArchRProj, 
    qc_metrix = c("Sample", "Clusters_qc", "log10(nFrags)", "DoubletEnrichment", "TSSEnrichment"), 
    embedding
    ){
    for (metrix in qc_metrix){
        # file_name <- 
        p <- plotEmbedding(
            ArchRProj = ArchRProj,
            colorBy = "cellColData",
            name = metrix,
            embedding = embedding)
        print(p)
    }
}


# 引数
# ArchRProject : proj
# ClustersName : ClustersAnno
# MatrixName : cisBPchromVar
# markerGenes : vectors of gene names
# motif : TRUE/FALSE
### TRUE  (for MotifMatrix)                : markerGenes -> markerMotifs, data -> z
### FALSE (for GeneScore, GeneIntegration) : markerGenes, data -> GeneScoreMatrix

# Object
### se : SummarizedExperiment obj.
### mat : Matrix obj.
### df_long : data.frame
### df_wide : data.frame
##### row     : gene names (for MotifMatrix, change names from motif_id to symbol)
##### columns : cell_types

# return
### df_wide : data.frame

getMatrixHeatmap <- function(
    ArchRProj,
    ClustersName,
    MatrixName,
    markerGenes,
    cellOrder,
    motif
) {
    se <- getMatrixFromProject(ArchRProj, MatrixName)
    col_meta <- getCellColData(ArchRProj)[[ClustersName]]
    row_meta <- se@elementMetadata$name

    if (motif == TRUE) {
        # get motif ids
        motifPositions <- getPositions(ArchRProj)
        markerGenes <- unlist(lapply(paste0("^", markerGenes, "_"), function(x) grep(x, names(motifPositions), value = TRUE)))
        
        ### -----
        mat <- se@assays@data$z
        colnames(mat) <- col_meta
        rownames(mat) <- row_meta

        df_long <- reshape2:::melt.matrix(mat)
        colnames(df_long) <- c("symbol", "cell_type", "score")

        df_wide <- df_long %>%
            group_by(cell_type, symbol) %>%
            summarize(mean_score = mean(score)) %>%
            spread(key = cell_type, value = mean_score) %>%
            as.data.frame()
        
        rownames(df_wide) <- df_wide$symbol
        df_wide <- df_wide %>% select(-symbol)

        out <- df_wide[markerGenes, cellOrder]
        out <- out  %>% 
            t() %>% 
            scale() %>% 
            t() %>%
            as.matrix()
        return(out)
        ### -----
        
    } else  {
        mat <- as.matrix(se@assays@data[[MatrixName]])
        colnames(mat) <- col_meta
        rownames(mat) <- row_meta

        df_long <- reshape2:::melt.matrix(mat)
        colnames(df_long) <- c("symbol", "cell_type", "score")

        df_wide <- df_long %>%
            group_by(cell_type, symbol) %>%
            summarize(mean_score = mean(score)) %>%
            spread(key = cell_type, value = mean_score) %>%
            as.data.frame()
        
        rownames(df_wide) <- df_wide$symbol
        df_wide <- df_wide %>% select(-symbol)

        out <- df_wide[markerGenes, cellOrder]
        out <- out  %>% 
            t() %>% 
            scale() %>% 
            t() %>%
            as.matrix()
        return(out)
    }
}