library(shiny)

source("/work/hello-testis/global.R")

#----------------------------------------------------------------------------------------------------
############## Server ##############
#----------------------------------------------------------------------------------------------------

shinyServer(function(input, output) {

    # eventReactive(input$goButton, {

    # Load ArchRProject


        output$UMAP <- renderPlot({
            # x <- faithful$waiting
            p <- plotEmbedding(
                proj, 
                colorBy = "cellColData", 
                embedding = "UMAP",
                name = "ClustersAnno")
                
            p1 <- plotEmbedding(
                ArchRProj = proj, 
                colorBy = "GeneScoreMatrix", 
                continuousSet = "horizonExtra",
                name = input$symbol, 
                embedding = "UMAP",
                imputeWeights = getImputeWeights(proj))

            p2 <- plotEmbedding(
                ArchRProj = proj, 
                colorBy = "GeneIntegrationMatrix", 
                name = input$symbol, 
                continuousSet = "coolwarm",
                embedding = "UMAP",
                imputeWeights = getImputeWeights(proj))

            ggAlignPlots(p, p1, p2, type = "h")

            # p2 <- plotBrowserTrack(
            #     ArchRProj = proj,
            #     groupBy = "ClustersAnno",
            #     useGroups = cellOrder,
            #     geneSymbol = c("Nanos2", "Stra8", "Sox9"),
            #     upstream = 250000,
            #     downstream = 250000,
            #     loops = getPeak2GeneLinks(proj))
        })

        ### Peaks to plot -----
        # features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"]
        # これをplotBrowserTrackに組み込む
        
        output$track <- renderPlot({
            # if (input$selected == NULL) {
                p <- plotBrowserTrack(
                    ArchRProj = proj, 
                    groupBy = "ClustersAnno", 
                    useGroups = cellOrder,
                    geneSymbol = input$symbol, 
                    upstream = 50000,
                    downstream = 50000,
                    loops = getPeak2GeneLinks(proj))
                # print(p)
            # } else {
            #     markersPeaks <- loadRDS("/path/to/markerPeaks.rds")
            #     p <- plotBrowserTrack(
            #         ArchRProj = proj, 
            #         groupBy = "ClustersAnno", 
            #         useGroups = cellOrder,
            #         geneSymbol = input$symbol, 
            #         features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[input$selected],
            #         upstream = 50000,
            #         downstream = 50000,
            #         loops = getPeak2GeneLinks(proj))
            #     print(p)
            # }

            grid::grid.newpage()
            grid::grid.draw(p[[input$symbol]])
        })

        output$etc <- renderPlot({

            motifs <- c(input$symbol)
            markerMotifs <- unlist(lapply(paste0("^", motifs, "_"), function(x) grep(x, names(motifPositions), value = TRUE)))

            p1 <- plotTrajectory(
                proj, 
                trajectory = "GermU", 
                colorBy = "GeneIntegrationMatrix", 
                name = input$symbol,
                imputeWeights = getImputeWeights(proj),
                addArrow = FALSE,
                continuousSet = "coolwarm")

            p2 <- plotTrajectory(
                proj, 
                trajectory = "GermU", 
                colorBy = "cisBPchromVar",
                name = paste0("z:", markerMotifs),
                imputeWeights = getImputeWeights(proj),
                addArrow = FALSE,
                continuousSet = "solarExtra")

            ggAlignPlots(p1[[2]], p2[[2]], type="h")
            # g1 <- ggplotGrob(p1[[2]])
            # g2 <- ggplotGrob(p2[[2]])
            # g <- rbind(g1, g2, size = "first")
            # g$widths = grid::unit.pmax(g1$widths, g2$widths)
        #     # plot(g)
        })

        output$footprint <- renderPlot({

            # if (containSoma) {
            #     cluster_name <- "ClustersAnno"
            # } else {
            #     cluster_name <- "ClustersAnno"  # Germ
            # }
            cluster_name <- "ClustersAnno"

            motifs <- c(input$symbol)
            markerMotifs <- unlist(lapply(paste0("^", motifs, "_"), function(x) grep(x, names(motifPositions), value = TRUE)))

            ### TF footprinting -----
            seFoot <- getFootprints(
                ArchRProj = proj, 
                positions = motifPositions[markerMotifs], 
                groupBy = cluster_name)
            
            p <- plotFootprints(
                seFoot = seFoot,
                ArchRProj = proj, 
                normMethod = "Subtract",
                # plotName = paste0("Footprints-Subtract-Bias-", motifs),
                addDOC = FALSE,
                smoothWindow = 5,
                baseSize = 12,
                height = 20,
                width = 6,
                plot = FALSE)
            # print(p$markerMotifs)
            # grid::grid.newpage()
            # grid::grid.draw(p$markerMotifs)
            gridExtra::grid.arrange(p[[markerMotifs]], heights=20, widths=16)
        })

})