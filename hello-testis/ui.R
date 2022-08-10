library(shiny)

genes <- readRDS("/work/hello-testis/genes.rds")

#----------------------------------------------------------------------------------------------------
############## UI ##############
#----------------------------------------------------------------------------------------------------

shinyUI(
    tagList(
        fluidPage(

        # App title ----
        titlePanel("Hello, testis",
            tags$head(tags$link(rel = "icon", type = "image/pmg", href = "logo.png"))
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
                
                selectInput(
                    "peaks",
                    label = "Select peaks set to plot (D)\n (default : all detected peaks)",
                    choices = c("all", cellOrder),
                    selected = "all" # default : all
                ),

                selectInput(
                    "motifs",
                    label = "Select motifs to plot \n (default : NONE)",
                    choices = names(motifPositions),
                    selected = "Dmrt1_251"
                ),

                checkboxInput(
                    "chip",
                    label = "Display ChIP-seq peaks (D)\n (default : FALSE)"
                ),

                checkboxInput(
                    "footprint",
                    label = "Display TF footprinting if the motif exists (G) (Plotting takes more time.)"
                ),

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