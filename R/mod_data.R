#' data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_data_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_data_ui(ns=",ns("foo"),")"))
  tagList(
    fluidRow(
      column(6,
             selectInput(inputId  = ns("data_source"),
                         label    = "Data source",
                         choices  = parse_data_source_choices(id),
                         selected = ""),
             sliderTextInput(inputId  = ns("pval"),
                             label    = "p-value thresh",
                             choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                             selected = 1.0,
                             grid     = TRUE),
             actionButton(inputId = ns("import"),
                          label   = "Import"),
      ),
      column(6,
             uiOutput(ns("source_selector"))
      ),
    )
  )
}


#' @title Parse UI data choices
#' @description
#' Takes the parent module UI `id` and then populates it's data module with only those
#' data sources that are relevant to that particular module. e.g. the GWAS module can
#' import from GWAS sources, but doesn't need access to the eQTL sources etc.
#' @return a character vector of choices to populate the `selectInput(inputId  = ns("data_source"))` UI element
#' @noRd
#'
parse_data_source_choices <- function(id) {

  # remove the module number (there may be multiple GWAS modules active - GWAS_0, GWAS_1 ... etc)
  # giving data modules GWAS_0-data, GWAS_1-data ... etc)
  # gotcha catch for wrong data id naming
  if(!grepl("([A-z])_\\d+-data$", id)) {
    warning("Data module: `", id, "` - please name all data modules as `[parent_id]-data` (replace square brackets and parent_id with module name)")
    stop()
  }

  # extract the parent ID
  parent_id_type <- sub("([A-z])_\\d+-data","\\1",id)

  # Return the data sources relevant to the module
  choices <- switch(parent_id_type,
                    "GWAS"  = c("Local","IEU Open GWAS","EBI GWAS Catalogue","EBI eQTL Catalogue"),
                    "eQTL"  = c("EBI eQTL Catalogue"),
                    "Coloc" = c("EBI eQTL Catalogue"))
}


#' data Server Functions
#' @importFrom httr GET stop_for_status content
#' @noRd
mod_data_server <- function(id, gene_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_data_server(ns=",ns("foo"),")"))

    #==========================================
    # imported dataset info tables and datasets
    #==========================================
    v <- reactiveValues(
      source = "",
      source_id = "",
      info   = NULL,
      data   = NULL,
      data2  = NULL
    )


    #==========================================
    # observe the data source select input
    #==========================================
    session$userData[[ns("data_source")]] <- observeEvent(input$data_source, {
      cli::cli_alert_info(paste0("mod_data_server::observeEvent(input$data_source)"))

      # deal with other types of data input; v$source check is to prevent repeated calls
      if (input$data_source=="EBI eQTL Catalogue" && v$source!="EBI eQTL Catalogue") {

        # query the available datasets from the EBI API
        v$info   <- get_ebi_eqtl_info()
        v$source <- "EBI eQTL Catalogue"

      } else if (input$data_source=="IEU Open GWAS" && v$source!="IEU Open GWAS") {

        # query the available datasets from the IEU API
        studies <- get_ieu_gwas_info()

        if(!is.data.frame(studies)) { # happens when things fail
          tries   <- 0
          while(!is.data.frame(studies) & tries < 5) {
            memoise::drop_cache(get_ieu_gwas_info)() # remove from cache
            studies <- get_ieu_gwas_info()           # try again
            tries <- tries + 1
          }
        }
        v$info   <- studies
        v$source <- "IEU Open GWAS"


      } else if (input$data_source=="EBI GWAS Catalogue" && v$source!="EBI GWAS Catalogue") {

        # query the available datasets from the EBI GWAS API
        v$info   <- get_ebi_gwas_info()
        v$source <- "EBI GWAS Catalogue"

      } else if (input$data_source=="Local" && v$source!="Local") {

        # query the available datasets from TargetExplorerData
        v$info   <- get_available_data_sources()
        v$source <- "Local"

      } else {

        # await file input as 4th option

      }

    })


    #==========================================
    # study information filtering module
    #==========================================
    output$source_selector <- renderUI({
      cli::cli_alert_info(paste0("mod_data_server::output$source_selector <- renderUI"))

      # TODO: fix the bug that removing and readding modules breaks auto source_selector UIoutput
      req(input$data_source)

      if(input$data_source %in% c("IEU Open GWAS","EBI GWAS Catalogue")) {

        selector <- selectizeInput(inputId = ns("filter_gwas"),
                                   label   = "GWAS ID",
                                   choices = NULL)

        # server side choices update much faster than client side (https://shiny.posit.co/r/articles/build/selectize/)
        updateSelectizeInput(session, 'filter_gwas', choices=v$info$id, server=TRUE)

        return(selector)

      } else if(input$data_source=="File input") {

        selector <- fileInput(inputId = ns("file"),
                              label   = "File",
                              buttonLabel = "file path...")

        return(selector)

      } else if(input$data_source == "EBI eQTL Catalogue") {

        selector <- datamods::select_group_ui(id = ns("filter"),
                                              params = list(
                                                catalogue    = list(inputId = "catalogue",   label = "Catalogue:"),
                                                dataset      = list(inputId = "dataset",     label = "Dataset:"),
                                                quant_method = list(inputId = "quant_method",label = "Quant. method:")
                                              ),
                                              vs_args = list(noOfDisplayValues = 1),
                                              inline = FALSE)

        filter <- datamods::select_group_server(id="filter",
                                                data = reactive(v$info),
                                                vars = reactive(c("catalogue","dataset","quant_method")))

        return(selector)

      } else if(input$data_source == "Local") {

        selector <- datamods::select_group_ui(id = ns("filter"),
                                              params = list(
                                                catalogue    = list(inputId = "catalogue",   label = "Catalogue:"),
                                                dataset      = list(inputId = "dataset",     label = "Dataset:")
                                              ),
                                              vs_args = list(noOfDisplayValues = 1),
                                              inline = FALSE)

        filter <- datamods::select_group_server(id="filter",
                                                data = reactive(v$info),
                                                vars = reactive(c("catalogue","dataset")))

        return(selector)

      }

    })


    #==========================================
    # observe the import button
    #==========================================
    session$userData[[ns("import")]] <- observeEvent(input$import, {

      # custom user file input
      if(input$data_source=="File input") {

        if(is.null(input$gwas_file)) {
          return()
        } else {
          # TODO
          # v$source_id <- input$filter_gwas
        }

      # OpenGWAS API
      } else if(input$data_source=="IEU Open GWAS") {

        v$source_id <- input$filter_gwas

        v$data <- import_ieu_gwas(study_id = input$filter_gwas,
                                  chr      = gene_module$chr,
                                  bp_start = gene_module$start - gene_module$flanks_kb*1000,
                                  bp_end   = gene_module$end + gene_module$flanks_kb*1000,
                                  build    = gene_module$build)

        # if fails, the memoise cache will have stored NULL, need to remove to enable a rerun next time
        if(is.null(v$data)) {
          memoise::drop_cache(import_ieu_gwas)(study_id = input$filter_gwas,
                                               chr      = gene_module$chr,
                                               bp_start = gene_module$start - gene_module$flanks_kb*1000,
                                               bp_end   = gene_module$end + gene_module$flanks_kb*1000,
                                               build    = gene_module$build)

          # if not null then also get some gene data
        } else {
          v$data2 <- read_gene_data(source = "gencode",
                                    chrom  = gene_module$chr,
                                    start  = gene_module$start - gene_module$flanks_kb*1000,
                                    end    = gene_module$end + gene_module$flanks_kb*1000)
        }

      # EBI GWAS catalogue API
      } else if(input$data_source=="EBI GWAS Catalogue") {

        v$source_id <- input$filter_gwas

        v$data <- import_ebi_gwas(study_id = input$filter_gwas,
                                  chr      = gene_module$chr,
                                  bp_start = gene_module$start - gene_module$flanks_kb*1000,
                                  bp_end   = gene_module$end + gene_module$flanks_kb*1000,
                                  build    = gene_module$build,
                                  pthresh  = input$pval)


      # EBI eQTL catalogue API
      } else if(input$data_source=="EBI eQTL Catalogue") {

        v$source_id <- paste0(input$`filter-dataset`, collapse = ", ")

        v$data <- import_ebi_eqtl(study_info_table = v$info,
                                  studies  = input$`filter-catalogue`,
                                  tissues  = input$`filter-dataset`,
                                  method   = input$`filter-quant_method`,
                                  chr      = gene_module$chr,
                                  bp_start = gene_module$start - gene_module$flanks_kb*1000,
                                  bp_end   = gene_module$end + gene_module$flanks_kb*1000,
                                  gene_id  = gene_module$id,
                                  pthresh  = input$pval,
                                  build    = gene_module$build)

        v$data2 <- read_gene_data(source = "gencode",
                                  chrom  = gene_module$chr,
                                  start  = gene_module$start - gene_module$flanks_kb*1000,
                                  end    = gene_module$end + gene_module$flanks_kb*1000)


      # Internal file
      } else if(input$data_source=="Local") {

        v$source_id <- paste0(input$`filter-dataset`, collapse = ", ")

        v$data <- read_internal_data(type    = input$`filter-catalogue`,
                                     dataset = input$`filter-dataset`,
                                     chrom   = gene_module$chr,
                                     start   = gene_module$start - gene_module$flanks_kb*1000,
                                     end     = gene_module$end + gene_module$flanks_kb*1000,
                                     pthresh = input$pval)

        v$data2 <- read_gene_data(source = "gencode",
                                  chrom  = gene_module$chr,
                                  start  = gene_module$start - gene_module$flanks_kb*1000,
                                  end    = gene_module$end + gene_module$flanks_kb*1000)

      }

    })


    # return the reactive data module values
    return(v)
  })
}

