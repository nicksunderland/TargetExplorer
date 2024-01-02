#' r_susier UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets shinyalert
#'
mod_r_susier_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_r_susier_ui(ns=",ns("foo"),")"))
  useShinyalert()
  tagList(
    fluidRow(
      column(6,
             selectInput(inputId  = ns("func"),
                         label    = "Function",
                         choices  = c("susie_rss"),
                         selected = "susie_rss")),
      column(4,
             actionButton(inputId = ns("run"),
                          label   = "Run")),
      column(1,
             actionButton(inputId = ns("reset"),
                          width   = "40px",
                          label   = "",
                          icon    = icon("rotate-left"))),
      tags$style(type='text/css', paste0("#",ns("run")," { width:100%; margin-top: 25px;}")),
      tags$style(type='text/css', paste0("#",ns("reset")," { width:100%; margin-top: 25px;}")),
    ),
    fluidRow(
      column(6,
             mod_reference_ui(id=ns("reference"))),
      column(3,
             prettyRadioButtons(inputId  = ns("source_type"),
                                label    = "Type",
                                choices  = c("quant","cc"),
                                selected = "quant",
                                inline   = FALSE)),
      column(3,
             numericInput(inputId  = ns("n_param"),
                          label    = "N",
                          value    = NA_real_,
                          step     = 1))
    ),
    fluidRow(
      column(3,
             numericInput(inputId = ns("sd_y"),
                          label = "sdY",
                          value = NA_real_,
                          step  = 0.1)),
      column(3,
             numericInput(inputId  = ns("s_param"),
                          label    = "S",
                          value    = NA_real_,
                          step     = 0.1)),
      column(3,
             numericInput(inputId  = ns("l_param"),
                          label    = "L",
                          value    = 5,
                          step     = 1)),
      column(3,
             numericInput(inputId  = ns("max_iter"),
                          label    = "Max iter",
                          value    = 100,
                          step     = 10))
    ),
    fluidRow(
      column(6,
             sliderTextInput(inputId  = ns("coverage"),
                             label    = "Coverage %",
                             choices  = c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0),
                             selected = 0.95,
                             grid     = TRUE)),
      column(6,
             sliderTextInput(inputId  = ns("maf_thresh"),
                             label    = "EAF thresh",
                             choices  = c(0.0,0.01,0.05,0.1,0.2),
                             selected = 0.0,
                             grid     = TRUE))
    )
  ) # tagList end
}


#' r_susier Server Functions
#' @noRd
mod_r_susier_server <- function(id, gene_module, data_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_r_susier_server(ns=",ns("foo"),")"))

    # R CMD checks
    # BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL
    run_flag <- reactiveVal(FALSE)

    #==========================================
    # Data, clump, and remove (sub)module servers for the GWAS module
    #==========================================
    reference_mod <- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=T)


    #==========================================
    # Run button
    #==========================================
    session$userData[[ns("run")]] <- observeEvent(input$run, {

      # check data
      if(is.null(data_module$data)) return(NULL)

      # if large amount of data, confirm run
      if(nrow(data_module$data) > 2000) {

        shinyalert(
          inputId = "shinyalert",
          title = "Processing time warning",
          text = paste0(nrow(data_module$data), " data points currently loaded, this may take some time. \n",
                        "Consider narrowing the genomic window and re-loading the data.\n",
                        "",
                        "Click OK to continue processing."),
          type = "warning",
          showCancelButton = TRUE,
          confirmButtonCol = "#AEDEF4",
          animation = FALSE,
          callbackR = function(ret) { if(ret) { run_flag(TRUE) } }
        )

      } else {

        run_flag(TRUE)

      }

    })


    #==========================================
    # Run function
    #==========================================
    session$userData[[ns("main")]] <- observeEvent(req(run_flag() == TRUE), {

      # reset run flag
      isolate(run_flag(FALSE))

      # reset if previously run (issues with factors being reset)
      if(any(c("group","index") %in% names(data_module$data))) {

        data_module$data$group <- NULL
        data_module$data$index <- NULL

      }

      # functions might fail
      tryCatch({

        shiny::withProgress(message = 'Finemapping data', value = 0, {

          shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))

          # finemap.abf
          if(input$func == "susie_rss") {

            # get the reference file
            if(grepl("1kGv3", reference_mod$ref_path)) {

              plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                           chrom    = gene_module$chr,
                                           from     = gene_module$start - gene_module$flanks_kb*1000,
                                           to       = gene_module$end + gene_module$flanks_kb*1000)

              dat_ld_obj <- genepi.utils::ld_matrix(dat       = data_module$data,
                                                    method    = "r",
                                                    plink2    = get_plink2_exe(),
                                                    plink_ref = plink_ref,
                                                    ukbb_ref  = NULL)

            } else if(grepl("UKB_LD", reference_mod$ref_path)) {

              # get the UKBB LD file path (downloads if not in cache)
              ukbb_ref_dt <- genepi.utils::download_ukbb_ld(chr           = gene_module$chr,
                                                            bp_start      = gene_module$start - gene_module$flanks_kb*1000,
                                                            bp_end        = gene_module$end + gene_module$flanks_kb*1000,
                                                            ukbb_ld_cache = file.path(reference_mod$ref_path, "cache"))

              # create the LD matrix for the region variants
              dat_ld_obj <- genepi.utils::ld_matrix(dat       = data_module$data,
                                                    method    = "r",
                                                    plink2    = get_plink2_exe(),
                                                    plink_ref = NULL,
                                                    ukbb_ref  = ukbb_ref_dt$root_file)
            }

            # store the r corr matrix
            data_module$ld_matrix_r  <- dat_ld_obj[["ld_mat"]]

            # reset the r2 matrix
            data_module$ld_matrix_r2 <- NULL

            # make and check the dataset
            D1 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
                                     type = input$source_type,
                                     sdY  = input$sd_y,
                                     S    = input$s_param,
                                     N    = input$n_param,
                                     ld   = dat_ld_obj[["ld_mat"]])

            # susie_rss
            results <- susieR::susie_rss(bhat  = D1$beta,
                                         shat  = D1$varbeta ^ 0.5,
                                         maf   = D1$MAF,
                                         R     = D1$LD,
                                         var_y = D1$sdY ^ 2,
                                         L     = input$l_param,
                                         n     = input$n_param,
                                         maf_thresh = input$maf_thresh,
                                         coverage = input$coverage,
                                         max_iter = input$max_iter)

            # update the data with grouping
            data_module$data <- calc_credible_set(result      = results,
                                                  dat_join_to = data_module$data,
                                                  coverage    = input$coverage,
                                                  result_type = "susie_rss")

            # store data for the kriging_rss plot
            data_module$kriging_rss         <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), D1$LD, n=D1$N)$conditional_dist
            data_module$kriging_rss$RSID    <- D1$snp
            data_module$kriging_rss$outlier <- ifelse(data_module$kriging_rss$logLR  > 2 &
                                                      abs(data_module$kriging_rss$z) > 2, TRUE, FALSE)


          }

          # warn if no groups found
          if(all(is.na(data_module$data$group))) {
            shiny::incProgress(3/4, detail = paste("Failed"))
            showNotification("Finemapping failed: \nno groups were found with the provided parameters, or susieR failed", type="error")
          } else {
            shiny::incProgress(3/4, detail = paste("Complete"))
          }

        }) # end withProgress

      },
      error=function(e) {

        showNotification(paste0(input$func, " failed - ", e), type="error", duration = 10)
        return(NULL)

      }) # end tryCatch

    }) # end Run button


    #==========================================
    # reset grouping
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
    # function / method choice and input control activation
    #==========================================
    observeEvent(list(input$func, input$source_type, data_module$data), {

      req(input$func)
      req(input$source_type)

      # common params
      if(input$source_type == "quant") {

        shinyjs::disable("s_param")
        shinyjs::enable("sd_y")

      } else if(input$source_type == "cc") {

        shinyjs::enable("s_param")
        shinyjs::disable("sd_y")

      }

      # sample size
      if(!is.null(data_module$data) && "N" %in% names(data_module$data) && is.na(input$n_param)) {

        sample_size <- max(data_module$data$N, na.rm=TRUE)
        updateNumericInput(inputId="n_param", value=sample_size)

      } else {

        sample_size <- input$n_param

      }

      # sdY
      if(!is.null(data_module$data) && all(c("SE","EAF") %in% names(data_module$data)) && !is.na(sample_size)) {

        if(input$source_type == "quant") {
          sdY_estimate <- sdY.est(data_module$data$SE ^ 2,
                                  data_module$data$EAF,
                                  sample_size)
          updateNumericInput(inputId="sd_y", value=sdY_estimate)
        } else {
          updateNumericInput(inputId="sd_y", value=NA_real_)
        }

      }

    }) # end controls setup


  })
}