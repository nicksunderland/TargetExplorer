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
  cli::cli_alert_info(paste0("initialising mod_grouping_ui(ns=",ns("foo"),")"))
  tagList(
    fluidRow(
      column(3, p(strong("Grouping:"), style="text-align:right; margin-top: 6px;")),
      column(9,
             selectInput(inputId  = ns("grouping"),
                         label    = NULL,
                         choices  = c("Off","plink::clump","r-coloc","r-susieR"),
                         selected = "Off"))
    ),
    uiOutput(ns("grouping_controls"))
  ) # tagList end
}


#' grouping Server Functions
#' @noRd
mod_grouping_server <- function(id, gene_module, data_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_grouping_server(ns=",ns("foo"),")"))

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL


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

      } else if(input$grouping == "plink::clump") {

        controls <- mod_clump_ui(id=ns("clump"))
        mod_clump_server(id="clump", gene_module, data_module)

      } else if(input$grouping == "r-coloc") {

        controls <- mod_r_coloc_ui(id=ns("r_coloc"))
        mod_r_coloc_server(id="r_coloc", gene_module, data_module)

      } else if(input$grouping == "r-susieR") {

        controls <- mod_r_susier_ui(id=ns("r_susier"))
        mod_r_susier_server(id="r_susier", gene_module, data_module)

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
