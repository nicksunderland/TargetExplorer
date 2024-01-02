#' reference UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_reference_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_reference_ui(ns=",ns("foo"),")"))
  tagList(
    uiOutput(ns("reference_ui"))
  )
}


#' reference Server Functions
#' @noRd
mod_reference_server <- function(id, gene_module, label="Reference", enabled=TRUE, include=NULL, exclude=NULL){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_reference_server(ns=",ns("foo"),")"))

    #==========================================
    # get the available genome references
    #==========================================
    ref_dir        <- get_data_directory(type="references")
    available_refs <- list.files(ref_dir, recursive = FALSE, full.names = TRUE)
    names(available_refs) <- basename(available_refs)
    if(!is.null(exclude)) {
      available_refs <- available_refs[!grepl(paste0(exclude,collapse="|"), available_refs)]
    }
    if(!is.null(include)) {
      available_refs <- available_refs[grepl(paste0(include,collapse="|"), available_refs)]
    }


    #==========================================
    # reactive values for this module
    #==========================================
    v <- reactiveValues(ref_path = NULL)


    #==========================================
    # reference selector UI component - # TODO - implement filtering of available reference depending on build type in `gene_module`
    #==========================================
    output$reference_ui <- renderUI({
      s <- selectInput(inputId = ns("reference"),
                       label   = label,
                       choices = c(names(available_refs)))
      if(enabled) {
        return(s)
      } else {
        return(shinyjs::disabled(s))
      }
    })


    #==========================================
    # observe the reference input selector
    #==========================================
    session$userData[[ns("reference")]] <- observeEvent(input$reference, {

      v$ref_path <- available_refs[[input$reference]]

    })


    #==========================================
    # Return the reference module reactive values
    #==========================================
    return(v)
  })
}



make_ref_subset <- function(ref_path, chrom, from, to) {

  # create file name for the reference subset
  cache_path <- file.path(ref_path, "cache")
  out <- file.path(cache_path, paste0(basename(ref_path),"_chr",chrom,"_",from,"_",to))

  # see if it already exists, if so return it
  if(paste0(out,".pgen") %in% list.files(cache_path, full.names=TRUE)) {
    return(out)
  } else {
    ref_path <- file.path(ref_path, basename(ref_path))
  }

  # else, create a subset of the reference file
  plink2 <- get_plink2_exe()

  # build the command
  cmd <- paste(plink2,
               "--pfile", ref_path, "vzs",
               "--chr", chrom,
               "--from-kb", floor(from/1000),
               "--to-kb", ceiling(to/1000),
               "--make-pgen", "vzs",
               "--out", out)

  # execute
  system(cmd)

  # confirm file
  if(!file.exists(paste0(out,".pgen"))) {
    stop("problems creating plink2 reference subset")
  } else {
    return(out)
  }
}



