#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # use shinyjs (for element activation/deactivation things)
    shinyjs::useShinyjs(),
    # Your application UI logic
    fluidPage(
      h1("TargetExplorer"),
      p("nicholas.sunderland@bristol.ac.uk"),
      hr(),
      fluidRow(
        column(10, mod_gene_ui(id="gene")),
        column(1, selectInput(inputId="module_type", label="Modules", choices=c("GWAS","eQTL","Coloc","MR"), selected = "GWAS")),
        column(1, actionButton(inputId="add_step", width = "40px", label = "", icon = icon("plus"))),
        tags$style(type='text/css', "#add_step { width:100%; margin-top: 25px;}")
      ),
      hr(),
      div(
        id = "module_placeholder"
      )
    )
  )
}


#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "TargetExplorer"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
