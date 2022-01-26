library(shiny)

genes <- readRDS("/work/hello-testis/genes.rds")

#----------------------------------------------------------------------------------------------------
############## UI ##############
#----------------------------------------------------------------------------------------------------

shinyUI(
    tagList(
        fluidPage(

        # App title ----
        titlePanel("Hello, testis"),

        # Sidebar layout with input and output definitions ----
        sidebarLayout(

            position = "left",

            # Sidebar panel for inputs ----
            sidebarPanel(

                helpText("From SSC(spermatogonial stem cells) to ES(elongating spermatids) during mouse spermatogenesis"),
                # img(src="/work/hello-testis/spermatogenesis.png"),

                # Input : Slider for the number of bins ----

                selectInput(
                    "symbol",
                    label = "Select a gene symbol to plot (B~F)",
                    choices = genes,
                    selected = "Stra8"
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
                    "chip",
                    label = "Display ChIP-seq peaks (D)\n (default : FALSE)"
                ),

                checkboxInput(
                    "containSoma",
                    label = "Containing Somatic cells (for TF footprinting) (G)\n (default : Only Germ cells)"
                ),
                
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

                h3("G. TF footprinting"))
        ),

        # Footer ----
        tags$footer(
            # tags$a("Maezawa lab", icon("github"))
            tags$a(href = "https://github.com/Hattyoriiiiiii", "Maezawa lab", icon("github"))
        ),

        includeCSS("style.css")
    )
))