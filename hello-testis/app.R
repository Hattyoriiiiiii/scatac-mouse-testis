library(shiny)
library(ArchR)
library(tidyverse)

set.seed(42)

print("Please wait...")
proj <- loadArchRProject(path = "Save-proj_traj/", showLogo = FALSE)
proj <- addImputeWeights(proj)

motifPositions <- getPositions(proj, name = "Motif_cisbp")
genes <- proj@geneAnnotation$genes$symbol %>% unique()
GermTrajectory <- c("Undiff_1", "Undiff_2", "Diff", "Lep", "Zyg_eP", "Pac_Dip", "RS_1", "RS_2", "RS_3", "ES")
cellOrder <- c("Undiff_1", "Undiff_2", "Diff", "Lep", "Zyg_eP", "Pac_Dip", "RS_1", "RS_2", "RS_3", "ES", "Sertoli", "Soma")
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersAnno", force=TRUE)


print(getAvailableMatrices(proj))

#----------------------------------------------------------------------------------------------------
############## UI ##############
#----------------------------------------------------------------------------------------------------

# Define UI for app that draws a histgram ----

ui <- tagList(
        fluidPage(

        # App title ----
        titlePanel("Hello, testis"),

        # Sidebar layout with input and output definitions ----
        sidebarLayout(

            position = "left",

            # Sidebar panel for inputs ----
            sidebarPanel(

                helpText("From Undifferented spermatogonia to elongating spermatids during mouse spermatogenesis"),

                # Input : Slider for the number of bins ----

                selectInput(
                    "symbol",
                    label = "Select a gene symbol to plot (B~F)",
                    choices = genes,
                    selected = "Yy1"
                    ),
                
                selectInput(
                    "peaks",
                    label = "Select peaks set to plot (D)\n (default : all detected peaks)",
                    choices = NULL,
                    selected = NULL # default : all
                ),

                selectInput(
                    "motifs",
                    label = "Select motifs to plot \n (default : NONE)",
                    choices = NULL,
                    selected = NULL
                ),

                checkboxInput(
                    "containSoma",
                    label = "Containing Somatic cells (for TF footprinting) (G)\n (default : Only Germ cells)"
                ),
                
                actionButton("goButton", "GO!!!")
            ),

            # Main panel for displaying outputs ----
            mainPanel(
                # img(src = "P4bBKfL.png", height = 140, width = 140)
                # Output: histgram ----
                h3("A-C. UMAP representations"),
                h4("A. Cell type, B. Gene score, C. Gene expression"),
                plotOutput(outputId = "UMAP"),
                h3("D. Genome browser representation"),
                plotOutput(outputId = "track"),
                h3("E-F. Trajectory"),
                h4("E. Gene expression along with pseudotime"),
                h4("F. Motif accessibility along with pseudotime"),
                plotOutput(outputId = "etc"),
                h3("G. TF footprinting")
            )
        ),

        tags$footer(
            # tags$a("Maezawa lab", icon("github"))
            tags$a(href = "https://github.com/Hattyoriiiiiii", "Maezawa lab", icon("github"))
        ),

        includeCSS("style.css")
    )
)

#----------------------------------------------------------------------------------------------------
############## Server ##############
#----------------------------------------------------------------------------------------------------

server <- function(input, output) {

    # observeEvent(input$goButton, {

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


        #     ### TF footprinting -----
            # seFoot <- getFootprints(
            #     ArchRProj = proj, 
            #     positions = motifPositions[markerMotifs], 
            #     groupBy = "ClustersAnno")
            
            # p <- plotFootprints(
            #     seFoot = seFoot,
            #     ArchRProj = proj, 
            #     normMethod = "Subtract",
            #     plotName = "Footprints-Subtract-Bias",
            #     addDOC = FALSE,
            #     smoothWindow = 5,
            #     baseSize = 12,
            #     height = 20,
            #     width = 6,
            #     plot = FALSE)
            # gridExtra::grid.arrange(g, p$markerMotifs, heights=c(20,2), widths=c(16,20))
        })
    # })
}

# p1 <- plotEmbedding()
# p2 <- plotEmbedding()
# gglignPlots(p1, p2, type = "h")

# query_genes <- c()
# p <- plotEmbedding()

shinyApp(ui = ui, server = server)