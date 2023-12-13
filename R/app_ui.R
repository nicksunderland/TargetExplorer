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
        column(6, mod_base_ui(id="base")),
        column(2, selectInput(inputId="module_type", label="Modules:", choices=c("GWAS","eQTL","Coloc","MR"))),
        column(1, actionButton(inputId="add_step", label="Add")),
        column(2, selectInput(inputId="active_modules", label="Current modules:", choices=c(""))),
        column(1, actionButton(inputId="remove_step", label="Remove")),
        tags$style(type='text/css', "#add_step { width:100%; margin-top: 25px;}"),
        tags$style(type='text/css', "#remove_step { width:100%; margin-top: 25px;}")
      ),
      hr(),
      div(
        id = "add_here"
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
