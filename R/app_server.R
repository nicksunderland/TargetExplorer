

#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # a list of the active modules
  active_modules <- reactiveValues(modules = list(
    "base" = mod_base_server(id="base")
  ))

  # observe the add step button
  observeEvent(input$add_step, {

    # get the module type to add and how many already exist
    add_what <- input$module_type
    how_many <- sum(grepl(input$module_type, names(active_modules$modules)))

    # generate a unique id; zero indexed
    this_id  <- paste0(input$module_type,"_",how_many)

    # add the correct module type
    if(add_what=="GWAS") {

      # create the module and add to the list of active modules
      active_modules$modules[[this_id]] <- mod_gwas_server(id           = this_id,
                                                           base_module  = active_modules$modules[['base']])

      # add the correct UI type
      insertUI(selector="#add_here", ui=mod_gwas_ui(id=this_id))

    } else if(add_what=="Coloc") {

      # create the module and add to the list of active modules
      active_modules$modules[[this_id]] <- mod_coloc_server(id           = this_id,
                                                           base_module  = active_modules$modules[['base']])

      # add the correct UI type
      insertUI(selector="#add_here", ui=mod_coloc_ui(id=this_id))

    } else if(add_what=="MR") {

    }

    # update the active_modules select box
    updateSelectInput(inputId="active_modules", choices=names(active_modules$modules)[names(active_modules$modules)!="base"])

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
    active_modules$modules[[unwanted_step_id]] <- NULL

    # update the active_modules select box
    updateSelectInput(inputId="active_modules", choices=names(active_modules$modules)[names(active_modules$modules)!="base"])
  })

}
