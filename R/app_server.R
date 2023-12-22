

#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # a list of the active modules
  app <- reactiveValues(modules = list(
    "gene"  = mod_gene_server(id="gene")
  ))

  # the types of module that can be added
  module_options <- list(
    'GWAS'  = list(svr='mod_gwas_server',  ui='mod_gwas_ui',  id_idx=0),
    'eQTL'  = list(svr='mod_eqtl_server',  ui='mod_eqtl_ui',  id_idx=0),
    'Coloc' = list(svr='mod_coloc_server', ui='mod_coloc_ui', id_idx=0),
    'MR'    = list(svr='mod_mr_server',    ui='mod_mr_ui',    id_idx=0)
  )

  # observe the add step button
  observeEvent(input$add_step, {
    cli::cli_alert_info("app_server::observeEvent(input$add_step)")

    # generate a unique id for the module e.g. 'GWAS_0'
    id <- paste0(input$module_type,"_",module_options[[input$module_type]][['id_idx']])

    # create the module server and add to the list of app modules
    app$modules[[id]] <- do.call(module_options[[input$module_type]][['svr']], list(id=id, app=app))

    # create the module ui
    ui <- do.call(module_options[[input$module_type]][['ui']], list(id=id))

    # add the module UI to the main UI
    insertUI(selector="#module_placeholder", ui=ui)

    # increment the index (I cant figure out how to remove all references to the id used to create
    # a module, so for now need to create total unique, even if the module has previously been removed)
    module_options[[input$module_type]][['id_idx']] <<- module_options[[input$module_type]][['id_idx']] + 1

  })



}
