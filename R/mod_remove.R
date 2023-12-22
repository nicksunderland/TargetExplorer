#' remove UI Function
#'
#' @description
#' A shiny sub-Module to remove the parent module it is placed within. The UI is
#' a simple bin icon button and the action is to remove the parent module's UI,
#' inputs, and observers.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_remove_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_remove_ui(ns=",ns("foo"),")"))
  tagList(
    actionButton(inputId=ns("remove_module"), width = "40px", label = "", icon = icon("trash-can"))
  )
}


#' remove Server Functions
#' @noRd
mod_remove_server <- function(id, app, parent_id, parent_inputs){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_remove_server(ns=",ns("foo"),") id=",id))

    #==========================================
    # observe the remove button
    #==========================================
    observeEvent(input$remove_module, {
      cli::cli_alert_info("app_server::observeEvent(input$remove_step)")

      # remove from the UI
      removeUI(
        selector = paste0("#",parent_id),
        session = session
      )

      # remove the shiny inputs for this id
      # this code is pulled from a custom function (see this blog:
      # https://appsilon.com/how-to-safely-remove-a-dynamic-shiny-module/
      lapply(paste0(parent_id,"-",names(parent_inputs)), function(i) { .subset2(parent_inputs, "impl")$.values$remove(i) })

      # destroy() then remove all the observers with a matching id - this is the bit I'm not sure if it actually works, as
      # deleting then creating (with the same id) a module doesnt call all the reactives at start up - there must be copies
      # of the initial id floating around somewhere that thinks it is already registered...?
      # **N.B**
      # for this to work all observers need to be added to the `session$userData` list, like so:
      # session$userData[[ns("data_source")]] <- observeEvent(.....)
      # *******
      module_observers <- grep(parent_id, names(session$userData), value=TRUE)
      lapply(module_observers, function(i) { session$userData[[i]]$destroy() })
      rm(list=ls(module_observers, envir=session$userData), envir=session$userData)

      # remove from the active_modules reactive list
      app$modules[[parent_id]] <- NULL

    })


  })
}

