#' source_select UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_source_select_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("source_select_ui"))
  )
}


#' source_select Server Functions
#' @noRd
mod_source_select_server <- function(id, app, source_type=c("GWAS","Coloc","MR"), label="Source 1"){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #==========================================
    # reactive values for this module
    #==========================================
    v <- reactiveValues(data = NULL,
                        source_id = "")


    #==========================================
    # source_select UI component
    #==========================================
    output$source_select_ui <- renderUI({

      # create the selectInput
      s <- selectInput(inputId = ns("source"),
                       label   = label,
                       choices = c("-"))

      # add an observer listening for changes to the app modules (has an init call to populate the initial UI element)
      session$userData[[ns("app-modules")]] <- observeEvent(app$modules, {

        # the names of all the active modules in the app `modules` reactive values list
        active_module_ids <- names(app$modules)

        # just the modules to consider importing from for this selectInput
        select_module_ids <- active_module_ids[grepl(paste0(source_type, collapse="|"), active_module_ids) & !active_module_ids %in% c("gene", "id")]

        # get the currently selected input source
        selected <- input$source

        # update the selectInput
        updateSelectInput(inputId = "source", choices = c("-", select_module_ids), selected = selected)

      })

      # return the UI element
      return(s)
    })


    #==========================================
    # observe the source select for changes to selected source
    #==========================================
    session$userData[[ns("source")]] <- observeEvent(input$source, {

      # if there is an input data source
      if(!is.null(input$source) && !is.null(app$modules[[input$source]])) {

        # add an observer to track its data to this module (which is just a alias really)
        session$userData[[ns("source_data")]] <- observeEvent(app$modules[[input$source]]$data, {

          v$data      <- app$modules[[input$source]]$data
          v$source_id <- app$modules[[input$source]]$source_id

        })

      }

    })


    #==========================================
    # Return the reference module reactive values
    #==========================================
    return(v)
  })
}

