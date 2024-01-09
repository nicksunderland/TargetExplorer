#' gene UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(2, textInput(inputId    = ns("id"),
                          label      = "Ensembl ID:",
                          value      = "ENSG00000112164",
                          placeholder= "ENSG00000112164",
                          width      = "100%")),
      column(1, actionButton(inputId = ns("id_search"),
                             label   = "",
                             width   = "40px",
                             icon    = icon("search"))),
      tags$style(type='text/css', paste0("#",ns("id_search")," { width:100%; margin-top: 25px;}")),
      column(1, selectInput(inputId  = ns("chr"),
                            label    = "Chr",
                            choices  = c(as.character(1:22),"X"),
                            selected = "6")),
      column(2, numericInput(inputId = ns("start"),
                             label   = "Start",
                             value   = 39016574,
                             min     = 1,
                             max     = 250000000,
                             step    = 1000)),
      column(2, numericInput(inputId = ns("end"),
                             label   = "End",
                             value   = 39055519,
                             min     = 1,
                             max     = 250000000,
                             step    = 1000)),
      column(2, numericInput(inputId = ns("flanks_kb"),
                             label   = "Flanks-kb",
                             value   = 100,
                             min     = 0,
                             max     = 1000)),
      column(2, prettyRadioButtons(
                             inputId = ns("build"),
                             label   = "Build",
                             choices = c("GRCh37","GRCh38"),
                             selected= "GRCh37",
                             inline  = TRUE)),
    )
  )
}

#' gene Server Functions
#' @importFrom httr GET stop_for_status content
#' @noRd
mod_gene_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    v <- reactiveValues(
      chr       = NULL,
      start     = NULL,
      end       = NULL,
      flanks_kb = NULL,
      id        = NULL,
      build     = NULL
    )

    # observe the gene description inputs
    observeEvent(input$chr,       { v$chr       <- input$chr })
    observeEvent(input$start,     { v$start     <- input$start })
    observeEvent(input$end,       { v$end       <- input$end })
    observeEvent(input$flanks_kb, { v$flanks_kb <- input$flanks_kb })
    observeEvent(input$id,        { v$id        <- input$id })
    observeEvent(input$build,     { v$build     <- input$build })

    # observe the search for gene info button
    observeEvent(input$id_search, {

      tryCatch({

        gene_dat <- get_ensembl_gene_info(input$id, input$build)

        # update the UI
        updateSelectInput( inputId = "chr",   selected = gene_dat$seq_region_name)
        updateNumericInput(inputId = "start", value    = gene_dat$start)
        updateNumericInput(inputId = "end",   value    = gene_dat$end)

      }, error=function(e) {

        updateTextInput(inputId="id", value="", placeholder=e$message)

      }) # end tryCatch

    })

    # return the reactive values
    return(v)
  })
}

