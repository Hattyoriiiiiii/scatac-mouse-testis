
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

