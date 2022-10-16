library(shiny)
library(shinybusy)

genes <- readRDS("/work/hello-testis/genes.rds")

#----------------------------------------------------------------------------------------------------
############## UI ##############
#----------------------------------------------------------------------------------------------------

shinyUI(
    tagList(
        fluidPage(

        add_busy_gif(
        src = "https://jeroen.github.io/images/banana.gif",
        height = 100, width = 100
        ),

        # App title ----
        titlePanel("Hello, testis",
            tags$head(tags$link(rel = "icon", type = "image/png", href = "logo.png"))
        ),

        # sidebarLayout(
        #     position = "right",
        #     helpText("The figure below is a good combination of Shogi pieces, just like scRNA and scATAC."),
        #     img(src="logo.png")
        # ),

        # Sidebar layout with input and output definitions ----
        sidebarLayout(

            position = "left",

            # Sidebar panel for inputs ----
            sidebarPanel(
                # img(src="logo.png"),
                helpText("From SSC(spermatogonial stem cells) to ES(elongating spermatids) during mouse spermatogenesis"),

                # Input : Slider for the number of bins ----

                selectInput(
                    "symbol",
                    label = "Select a gene symbol to plot (B~F)",
                    choices = genes,
                    selected = "Sox4"
                    ),

                sliderInput("distance", label = "Distance (kb):", min = -500, max = 500, value = c(-50, 50)),

                selectInput(
                    "peaks",
                    label = "Select peaks set to plot (D)\n (default : all detected peaks)",
                    choices = c("all", cellOrder),
                    selected = "all" # default : all
                ),

                selectInput(
                    "tfbs_tf",
                    label = "Select TFBS to plot \n (default : Dmrt1)",
                    choices = c(TFBSs),
                    selected = "Dmrt1_MA1603.1"
                ),

                selectInput(
                    "tfbs_celltype",
                    label = "Select cell-type for TFBS to plot \n (default : SSC_1)",
                    choices = cellOrder,
                    selected = "SSC_1"
                ),

                # selectInput(
                #     "motifs",
                #     label = "Select motifs to plot \n (default : Dmrt1)",
                #     choices = names(motifPositions),
                #     selected = "Dmrt1_251"
                # ),

                selectInput(
                    "chip_tf",
                    label = "Display ChIP-seq peaks (D)\n (default : DMRT1)",
                    choices = c(chip),
                    selected = "DMRT1"
                ),

                # checkboxInput(
                #     "footprint",
                #     label = "Display TF footprinting if the motif exists (G) (Plotting takes more time.)"
                # ),

                # checkboxInput(
                #     "containSoma",
                #     label = "Containing Somatic cells (for TF footprinting) (G)\n (default : Only Germ cells)"
                # ),
                
                actionButton("goButton", "GO!!!")),


            # Main panel for displaying outputs ----
            mainPanel(
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

                h3("G. TF footprinting"),
                plotOutput(outputId = "footprint")
                ),

        ),

        # Footer ----
        tags$footer(
            tags$a(href = "https://github.com/Hattyoriiiiiii", "Maezawa lab", icon("github"))
        ),

        includeCSS("www/style.css")
    )
))


# about_panel <- tabPanel(

#   titlePanel(h5("About")),
#   # Various tabs.
#   tabsetPanel(
#     # General info.
#     tabPanel(
#       "Overview",
#       tags$h3("Scope"),
#       tags$p(HTML("ShinyArchR.UiO is a user-friendly, integrative open-source Shiny-based web app using R programming for visualization of massive single-cell chromatin accessibility data (scATAC-seq) based on <a href=\"https://www.archrproject.com\" target=\"_blank\">ArchR</a> (Corces et al., 2021).")),
#       tags$h3("Approach"),
      
#       tags$p(HTML(" The ArchR objects saved in folders along with HDF5 formatted Arrow files are used for input in ShinyArchR.UiO.")),
#       tags$h5("Data Visualization of ShinyArchR.UiO:"),
#       tags$ul(
#         tags$li(HTML("scATAC-seq clusters, unconstrained and constrained clusters on integrated reduced dimensions UMAP from ArchR objects")),
#         tags$li(HTML("Peaks using plot browser tracks on clusters on scATAC-seq modality")),
#         tags$li(HTML("Peaks2Genelinks tracks on single-cell RNA sequencing (scRNA-seq) integrated data with scATAC-seq using plot browser tracks. The co-accessibility among the genes can be visualized in the bottom panel")),
#         tags$li(HTML("Heatmaps of pseudo time trajectory")),
#         tags$li(HTML("Heatmaps of top 50 markers Peak2Genelinks in scATAC-seq and scRNA-seq")),

#         )
#     ),
    

#     # About us .
#     tabPanel(
#       "About us ",
#       tags$h3("Contributions and Citation info"),
#       tags$p(HTML("ShinyArchR.UiO software is developed by <a href=\"https://github.com/Hattyoriiiiiii\" target=\"_blank\">Hattyoriiiiiii</a> at <a href=\"https://www.uio.no\" target=\"_blank\">Tokyo University of Science</a>, as an open-source project mainly under the GPL license version 3 (see source code for details).")),
#       tags$p(HTML("If ShinyArchR.UiO in any way help you in visualizing and sharing your research work such that it warrants a citation, please cite the ShinyArchR.UiO preprint in BioRXiv or the final publication.")),
#     )

#   )

# )