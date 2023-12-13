#' base UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_base_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3, selectInput(inputId  = ns("gene_chr"),
                            label    = "Chrom:",
                            choices  = c(as.character(1:22),"X"),
                            selected = "6")),
      column(3, numericInput(inputId = ns("gene_start"),
                             label   = "Start:",
                             value   = 39016574,
                             min     = 1,
                             max     = 250000000,
                             step    = 1000)),
      column(3, numericInput(inputId = ns("gene_end"),
                             label   = "End:",
                             value   = 39055519,
                             min     = 1,
                             max     = 250000000,
                             step    = 1000)),
      column(3, numericInput(inputId = ns("gene_flanks_kb"),
                             label   = "Flanks (kb):",
                             value   = 250,
                             min     = 0,
                             max     = 1000))
    )
  )
}

#' base Server Functions
#'
#' @noRd
mod_base_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    base_vals <- reactiveValues(
      gene_chr       = NULL,
      gene_start     = NULL,
      gene_end       = NULL,
      gene_flanks_kb = NULL,
      data           = NULL
    )

    observeEvent(input$gene_chr,       { base_vals$gene_chr <- input$gene_chr })
    observeEvent(input$gene_start,     { base_vals$gene_start <- input$gene_start })
    observeEvent(input$gene_end,       { base_vals$gene_end <- input$gene_end })
    observeEvent(input$gene_flanks_kb, { base_vals$gene_flanks_kb <- input$gene_flanks_kb })

    return(base_vals)
  })
}

## To be copied in the UI
# mod_base_ui("base_1")

## To be copied in the server
# mod_base_server("base_1")
