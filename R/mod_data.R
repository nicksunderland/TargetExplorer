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
  tagList(
    fluidRow(
      column(6, selectInput(inputId  = ns("data_source"),
                            label    = "Data source",
                            choices  = c("Local","IEU Open GWAS","EBI eQTL Catalogue","File input"),
                            selected = "Local"),
             fileInput(inputId = ns("file"),
                       label   = "File",
                       buttonLabel = "file path..."),
             actionButton(inputId = ns("import"),
                          label   = "Import"),
      ),
      column(6,
             datamods::select_group_ui(id = ns("filter"),
                                       params = list(
                                         catalogue = list(inputId = "catalogue", label = "Catalogue:"),
                                         dataset   = list(inputId = "dataset",   label = "Dataset:")
                                       ),
                                       inline = FALSE),
             sliderTextInput(inputId  = ns("pval"),
                             label    = "p-value thresh",
                             choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                             selected = 1.0,
                             grid     = TRUE)
      ),
    )
  )
}

#' data Server Functions
#' @importFrom httr GET stop_for_status content
#' @noRd
mod_data_server <- function(id, gene){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    #==========================================
    # imported dataset info tables and datasets
    #==========================================
    v <- reactiveValues(
      source = "",
      info   = NULL,
      data   = NULL,
      data2  = NULL
    )


    #==========================================
    # study information filtering module
    #==========================================
    filter <- datamods::select_group_server(id="filter",
                                            data = reactive(v$info),
                                            vars = reactive(c("catalogue","dataset")))


    #==========================================
    # observe the data source select input
    #==========================================
    observeEvent(input$data_source, {

      # en/disable file input button
      if(input$data_source!="file input") {
        shinyjs::disable(id="file")
      } else {
        shinyjs::enable(id="file")
      }

      # deal with other types of data input
      if (input$data_source=="EBI eQTL Catalogue" && v$source!="EBI eQTL Catalogue") {

        # query the available datasets from the API
        v$info   <- get_ebi_study_info()
        v$source <- "EBI eQTL Catalogue"

      } else if (input$data_source=="IEU Open GWAS" && v$source!="EBI eQTL Catalogue" ) {

        # TODO

      } else if (input$data_source=="Local") {

        # query the available datasets from TargetExplorerData
        v$info   <- get_available_data_sources()
        v$source <- "Local"

      }

    })


    #==========================================
    # observe the import button
    #==========================================
    observeEvent(input$import, {

      # custom user file input
      if(input$data_source=="File input") {

        if(is.null(input$gwas_file)) {
          return()
        } else {
          # TODO
        }

      # OpenGWAS API
      } else if(input$data_source=="IEU Open GWAS") {

        # TODO

      # EBI eQTL catalogue API
      } else if(input$data_source=="EBI eQTL Catalogue") {

        v$data <- import_ebi_eqtl(study_info_table = v$info,
                                  studies  = input$`filter-catalogue`,
                                  tissues  = input$`filter-dataset`,
                                  chr      = gene$chr,
                                  bp_start = gene$start - gene$flanks_kb*1000,
                                  bp_end   = gene$end + gene$flanks_kb*1000,
                                  gene_id  = gene$id,
                                  nlog10p  = NULL, #-log10(input$pval),
                                  build    = gene$build)


      # internal file input
      } else if(input$data_source=="Local") {

        shiny::withProgress(message = 'Reading package data', value = 0, {
          shiny::incProgress(1/4, detail = paste("Dataset", input$`filter-dataset`, "..."))

          v$data <- read_internal_data(type    = input$`filter-catalogue`,
                                       dataset = input$`filter-dataset`,
                                       chrom   = gene$chr,
                                       start   = gene$start - gene$flanks_kb*1000,
                                       end     = gene$end + gene$flanks_kb*1000)

          v$data2 <- read_gene_data(source = "gencode",
                                    chrom  = gene$chr,
                                    start  = gene$start - gene$flanks_kb*1000,
                                    end    = gene$end + gene$flanks_kb*1000)

          shiny::incProgress(3/4, detail = paste("Complete"))
        })
      }

    })



    # return the reactive data module values
    return(v)
  })
}

