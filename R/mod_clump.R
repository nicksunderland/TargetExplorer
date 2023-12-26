#' clump UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_clump_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_clump_ui(ns=",ns("foo"),")"))
  tagList(
    fluidRow(
      column(6,
             mod_reference_ui(id=ns("reference"))),
      column(4,
             actionButton(inputId = ns("clump"),
                          label   = "Clump data")),
      column(1,
             actionButton(inputId = ns("reset"),
                          width   = "40px",
                          label   = "",
                          icon    = icon("rotate-left"))),
      tags$style(type='text/css', paste0("#",ns("clump")," { width:100%; margin-top: 25px;}")),
      tags$style(type='text/css', paste0("#",ns("reset")," { width:100%; margin-top: 25px;}")),
    ),
    uiOutput(ns("clump_params"))
  ) # tagList end
}


#' clump Server Functions
#' @noRd
mod_clump_server <- function(id, gene_module, data_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_gwas_server(ns=",ns("foo"),")"))

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL


    #==========================================
    # Data, clump, and remove (sub)module servers for the GWAS module
    #==========================================
    reference_mod <- mod_reference_server(id="reference", label="Clumping", gene_module=app$modules$gene)


    #==========================================
    # reference selector / controls
    #==========================================
    output$clump_params <- renderUI({

      if(!is.null(reference_mod$ref_path)) {

        shinyjs::enable("clump")
        shinyjs::enable("reset")

        ui_fluid_row <- fluidRow(
          column(6,
                 sliderTextInput(inputId  = ns("clump_p1"),
                                 label    = "p1",
                                 choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                                 selected = 5e-8,
                                 grid     = TRUE),
                 sliderTextInput(inputId  = ns("clump_r2"),
                                 label    = "r2",
                                 choices  = c(0.0001,0.001,0.01,seq(0.1,0.9,0.1),0.9999),
                                 selected = 0.001,
                                 grid     = TRUE)),
          column(6,
                 sliderTextInput(inputId  = ns("clump_p2"),
                                 label    = "p2",
                                 choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                                 selected = 1.0,
                                 grid     = TRUE),
                 sliderTextInput(inputId  = ns("clump_kb"),
                                 label    = "kb",
                                 choices  = c(1,seq(0,1000,50)[-1]),
                                 selected = 250,
                                 grid     = TRUE))
        )

      } else {

        shinyjs::disable("clump")
        shinyjs::disable("reset")

        # remove if previously run
        if("clump" %in% names(data_module$data)) {
          data_module$data$clump <- NULL
          data_module$data$index <- NULL
        }

        ui_fluid_row <- fluidRow()

      }

      return(ui_fluid_row)
    })


    #==========================================
    # reset clumping
    #==========================================
    session$userData[[ns("reset")]] <- observeEvent(input$reset, {

      if(all(c("clump","index") %in% names(data_module$data))) {
        data_module$data[, c("clump","index") := NULL]
        tmp <- data.table::copy(data_module$data)
        data_module$data <- NULL
        data_module$data <- tmp
      }

    })


    #==========================================
    # Clump button
    #==========================================
    session$userData[[ns("clump")]] <- observeEvent(input$clump, {

      # check data
      if(is.null(data_module$data)) return(NULL)

      # reset if previously run (issues with factors being reset)
      if(any(c("clump","index") %in% names(data_module$data))) {

        data_module$data$clump <- NULL
        data_module$data$index <- NULL

      }

      # get the reference file
      plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                   chrom    = gene_module$chr,
                                   from     = gene_module$start - gene_module$flanks_kb*1000,
                                   to       = gene_module$end + gene_module$flanks_kb*1000)


      shiny::withProgress(message = 'Clumping data', value = 0, {

        shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))

        # run clumping
        data_module$data <- genepi.utils::clump(gwas      = data_module$data,
                                                p1        = input$clump_p1,
                                                p2        = input$clump_p2,
                                                r2        = input$clump_r2,
                                                kb        = input$clump_kb,
                                                plink2    = get_plink2_exe(),
                                                plink_ref = plink_ref) #|> as.data.frame()

        # warn if no clumps found
        if(all(is.na(data_module$data$clump))) {
          shiny::incProgress(3/4, detail = paste("Failed"))
          showNotification("Clumping failed: \nno clumps were found with the provided parameters, or plink failed", type="error")
        } else {
          shiny::incProgress(3/4, detail = paste("Complete"))
        }

      })

    })

  })
}
