#' grouping UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_grouping_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3, p(strong("Grouping:"), style="text-align:right; margin-top: 6px;")),
      column(6,
             selectInput(inputId  = ns("grouping"),
                         label    = NULL,
                         choices  = c("Off","simple","plink::clump","r-coloc","r-susieR"),
                         selected = "Off"))
    ),
    uiOutput(ns("grouping_controls"))
  ) # tagList end
}


#' grouping Server Functions
#' @noRd
mod_grouping_server <- function(id, gene_module, data_module, parent_ui){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #==========================================
    # Reactive values
    #==========================================
    v <- reactiveValues(grouping = NULL)


    #==========================================
    # grouping selector
    #==========================================
    observeEvent(input$grouping, {

      v$grouping <- input$grouping

    })


    #==========================================
    # selector / controls
    #==========================================
    output$grouping_controls <- renderUI({

      if(input$grouping == "Off") {

        controls <- fluidRow()

        # reset if previously run (issues with factors being reset)
        if(!is.null(data_module$data) && any(c("group","index") %in% names(data_module$data))) {

          data_module$data$group <- NULL
          data_module$data$index <- NULL

        }

      } else if(input$grouping == "simple") {

        controls <- fluidRow(
          column(6, selectInput(inputId=ns("selection"), label="Select", choices=c("Off","All index","All group 1"), selected="Off"))
        )

        session$userData[[ns("data_update")]] <- observeEvent(data_module$data, {

          if("TISSUE" %in% names(data_module$data)) {
            current <- input$selection
            updateSelectInput(inputId = "selection", selected = current, choices =c("Off","All index","All group 1","TISSUE") )
          }

        })

        session$userData[[ns("selection")]] <- observeEvent(input$selection, {

          req(input$selection)
          req(data_module$data)

          if(input$selection=="Off"){
            if(any(c("group","index") %in% names(data_module$data))){
              data_module$data$group <- NULL
              data_module$data$index <- NULL
            }
          } else if(input$selection=="All index") {
            tmp <- data_module$data[,index:=TRUE]
            data_module$data <- NULL
            data_module$data <- tmp
          } else if(input$selection=="All group 1") {
            tmp <- data_module$data[,group:=factor(1)]
            data_module$data <- NULL
            data_module$data <- tmp
          } else{

            tmp <- data_module$data[,group:=factor(get(input$selection))]
            data_module$data <- NULL
            data_module$data <- tmp
          }

        })



      } else if(input$grouping == "plink::clump") {

        controls <- mod_clump_ui(id=ns("clump"))
        mod_clump_server(id="clump", gene_module=gene_module, data_module=data_module)

      } else if(input$grouping == "r-coloc") {

        controls <- mod_r_coloc_ui(id=ns("r_coloc"))
        mod_r_coloc_server(id              = "r_coloc",
                           gene_module     = gene_module,
                           source_1_module = data_module,
                           parent_ui       = parent_ui,
                           functions       = c("finemap.abf","finemap.signals"))

      } else if(input$grouping == "r-susieR") {

        controls <- mod_r_susier_ui(id=ns("r_susier"))
        mod_r_susier_server(id="r_susier", gene_module, data_module, parent_ui=parent_ui)

      }

      # return the controls UI elements
      return(controls)
    })


    #==========================================
    # Return the the active values for this module
    # just indicates whether it is in use/active or
    # or not.
    #==========================================
    return(v)
  })
}
