

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



  # observe the add step button
  observeEvent(input$add_step, {

    # get the module type to add and how many already exist
    add_what <- input$module_type
    how_many <- sum(grepl(input$module_type, names(app$modules)))

    # generate a unique id; zero indexed
    this_id  <- paste0(input$module_type,"_",how_many)

    # add the correct module type
    if(add_what=="GWAS") {

      # create the module and add to the list of active modules
      app$modules[[this_id]] <- mod_gwas_server(id = this_id, app = app)

      # add the correct UI type
      insertUI(selector="#add_here", ui=mod_gwas_ui(id=this_id))

    } else if(add_what=="Coloc") {

      # create the module and add to the list of active modules
      app$modules[[this_id]] <- mod_coloc_server(id = this_id, app = app)

      # add the correct UI type
      insertUI(selector="#add_here", ui=mod_coloc_ui(id=this_id))

    } else if(add_what=="MR") {

      # create the module and add to the list of active modules
      app$modules[[this_id]] <- mod_mr_server(id = this_id, app = app)

      # add the correct UI type
      insertUI(selector="#add_here", ui=mod_mr_ui(id=this_id))

    }

    # update the active_modules select box
    updateSelectInput(inputId = "active_modules",
                      choices = names(app$modules)[!names(app$modules) %in% c("gene")])

    names(app$modules) <- names(app$modules)

  })


  # observe the remove step button
  observeEvent(input$remove_step, {

    # get the id to remove
    unwanted_step_id <- input$active_modules

    # remove from the UI
    removeUI(
      selector = paste0("#",unwanted_step_id),
      session = session
    )

    # remove the shiny inputs
    remove_shiny_inputs(
      id = unwanted_step_id,
      .input = input
    )

    # may need to look into also removing the observer events
    # https://appsilon.com/how-to-safely-remove-a-dynamic-shiny-module/

    # remove from the active_modules reactive list
    app$modules[[unwanted_step_id]] <- NULL

    # update the active_modules select box
    updateSelectInput(inputId="active_modules",
                      choices=names(app$modules)[!names(app$modules) %in% c("gene")])
  })

}
