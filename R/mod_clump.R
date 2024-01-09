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
    fluidRow(
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
  ) # tagList end
}


#' clump Server Functions
#' @noRd
mod_clump_server <- function(id, gene_module, data_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL


    #==========================================
    # Data, clump, and remove (sub)module servers for the GWAS module
    #==========================================
    reference_mod <- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, exclude=c("UKB_LD"))


    #==========================================
    # reset clumping
    #==========================================
    session$userData[[ns("reset")]] <- observeEvent(input$reset, {

      if(all(c("group","index") %in% names(data_module$data))) {
        data_module$data[, c("group","index") := NULL]
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
      if(any(c("group","index") %in% names(data_module$data))) {

        data_module$data$group <- NULL
        data_module$data$index <- NULL

      }

      shiny::withProgress(message = 'Clumping data', value = 0, {

        shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))

        # get the reference file
        plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                     chrom    = gene_module$chr,
                                     from     = gene_module$start - gene_module$flanks_kb*1000,
                                     to       = gene_module$end + gene_module$flanks_kb*1000)

        # run clumping
        data_module$data <- genepi.utils::clump(gwas      = data_module$data,
                                                p1        = input$clump_p1,
                                                p2        = input$clump_p2,
                                                r2        = input$clump_r2,
                                                kb        = input$clump_kb,
                                                plink2    = get_plink2_exe(),
                                                plink_ref = plink_ref) #|> as.data.frame()

        # rename to standard index/group
        data.table::setnames(data_module$data, c("index","clump"), c("index","group"))

        # store the LD structure of the region
        data_module$ld_matrix_r2 <- genepi.utils::ld_matrix(dat       = data_module$data,
                                                            method    = "r2",
                                                            plink2    = get_plink2_exe(),
                                                            plink_ref = plink_ref)[["ld_mat"]]
        # store the LD structure of the region
        data_module$ld_matrix_r  <- genepi.utils::ld_matrix(dat       = data_module$data,
                                                            method    = "r",
                                                            plink2    = get_plink2_exe(),
                                                            plink_ref = plink_ref)[["ld_mat"]]

        # warn if no clumps found
        if(all(is.na(data_module$data$group))) {
          shiny::incProgress(3/4, detail = paste("Failed"))
          showNotification("Clumping failed: \nno clumps were found with the provided parameters, or plink failed", type="error")
        } else {
          shiny::incProgress(3/4, detail = paste("Complete"))
        }

      })

    })

  })
}
