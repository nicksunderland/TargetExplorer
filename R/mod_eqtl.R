#' eQTL UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_eqtl_ui <- function(id){
  ns <- NS(id)
  # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
  div(
    id = id,
    sidebarLayout(position = "right",
                  sidebarPanel(p(strong(paste0("Controls [",id,"]"))),
                               width = 3,
                               hr(),
                               mod_data_ui(id=paste0(id,"-data")),
                               hr(),
                               fluidRow(
                                 column(6,
                                        selectInput(inputId = ns("source_2"),
                                                    label   = "Upstream data",
                                                    choices = c(""))
                                 ),
                                 column(6,
                                        selectInput(inputId = ns("source_2_filter"),
                                                    label   = "Filter",
                                                    choices = c("All"),
                                                    selected = "All")
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 column(6,
                                        selectInput(inputId = ns("tissues"),
                                                    label   = "Tissues",
                                                    choices = c("All")),
                                        prettyRadioButtons(inputId = ns("beta_dir"),
                                                           label   = "BETA dir",
                                                           choices = c("Any","Concordant","Disconcordant"),
                                                           selected = "Any",
                                                           inline   = FALSE)
                                 ),
                                 column(6,
                                        selectInput(inputId = ns("tissue_combine_method"),
                                                    label   = "Combine by",
                                                    choices = c("Maximum BETA","Median BETA", "Lowest P-value", "Largest sample size"),
                                                    selected = "Largest sample size"),
                                        prettyRadioButtons(inputId = ns("output_variants"),
                                                           label   = "Output variants",
                                                           choices = c("All","Criteria"),
                                                           selected = "Criteria",
                                                           inline   = TRUE)
                                 )
                               ),

                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            fluidRow(
                              column(1, p(strong("eQTL:"))),
                              column(4, selectInput(inputId = ns("view"),
                                                    label   = NULL,
                                                    choices = c("Tissue overview", "Variant select"),
                                                    selected = "Tissue overview")),
                            ),
                            fluidRow(
                              column(6,
                                     plotOutput(outputId = ns("mr_plot"),
                                                height   = "500px",
                                                brush    = ns("mr_plot_brush")),
                                     tableOutput(outputId = ns("mr_result"))
                              )
                            ),
                  ) # main panel end
    ), # sidebar layout end
    hr()
  ) # div end
}


#' eQTL Server Functions
#' @import ggplot2 ggrepel
#' @noRd
mod_eqtl_server <- function(id, app){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # R CMD checks
    #

    #==========================================
    # Data module server for the eQTL module
    #==========================================
    data_mod <- mod_data_server(id="data", gene_module=app$modules$gene)


    #==========================================
    # eQTL plot
    #==========================================




    #==========================================
    # Return the data module to the reactive
    # values in the main app server such that
    # other modules can access the data to use
    # in their own processes.
    #==========================================
    return(data_mod)
  })
}
