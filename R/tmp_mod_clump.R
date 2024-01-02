#' #' clump UI Function
#' #'
#' #' @description A shiny Module.
#' #'
#' #' @param id,input,output,session Internal parameters for {shiny}.
#' #'
#' #' @noRd
#' #' @import shiny shinyWidgets
#' #'
#' mod_clump_ui <- function(id){
#'   ns <- NS(id)
#'   cli::cli_alert_info(paste0("initialising mod_clump_ui(ns=",ns("foo"),")"))
#'   tagList(
#'     fluidRow(
#'       column(3, p(strong("Grouping:"))),
#'       column(9,
#'              selectInput(inputId  = ns("grouping"),
#'                          label    = NULL,
#'                          choices  = c("Off","Clump","Bayesian"),
#'                          selected = "Off"))
#'     ),
#'     uiOutput(ns("grouping_controls"))
#'   ) # tagList end
#' }
#'
#'
#' #' clump Server Functions
#' #' @noRd
#' mod_clump_server <- function(id, gene_module, data_module){
#'   moduleServer( id, function(input, output, session){
#'     ns <- session$ns
#'     cli::cli_alert_info(paste0("initialising mod_gwas_server(ns=",ns("foo"),")"))
#'
#'     # R CMD checks
#'     BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL
#'
#'
#'     #==========================================
#'     # Reactive values
#'     #==========================================
#'     v <- reactiveValues(grouping = NULL)
#'
#'
#'     #==========================================
#'     # Data, clump, and remove (sub)module servers for the GWAS module
#'     #==========================================
#'     reference_mod <- mod_reference_server(id="reference", label="Reference", gene_module=app$modules$gene)
#'
#'
#'     #==========================================
#'     # grouping selector
#'     #==========================================
#'     observeEvent(input$grouping, {
#'
#'       v$grouping <- input$grouping
#'
#'     })
#'
#'
#'     #==========================================
#'     # selector / controls
#'     #==========================================
#'     output$grouping_controls <- renderUI({
#'
#'       if(input$grouping == "Off") {
#'
#'         controls <- fluidRow()
#'
#'         # reset if previously run (issues with factors being reset)
#'         if(!is.null(data_module$data) && any(c("group","index") %in% names(data_module$data))) {
#'
#'           data_module$data$group <- NULL
#'           data_module$data$index <- NULL
#'
#'         }
#'
#'       } else if(input$grouping == "Clump") {
#'
#'         controls <- tagList(
#'           fluidRow(
#'             column(6,
#'                    mod_reference_ui(id=ns("reference"))),
#'             column(4,
#'                    actionButton(inputId = ns("clump"),
#'                                 label   = "Clump data")),
#'             column(1,
#'                    actionButton(inputId = ns("reset"),
#'                                 width   = "40px",
#'                                 label   = "",
#'                                 icon    = icon("rotate-left"))),
#'             tags$style(type='text/css', paste0("#",ns("clump")," { width:100%; margin-top: 25px;}")),
#'             tags$style(type='text/css', paste0("#",ns("reset")," { width:100%; margin-top: 25px;}")),
#'           ),
#'           fluidRow(
#'             column(6,
#'                    sliderTextInput(inputId  = ns("clump_p1"),
#'                                    label    = "p1",
#'                                    choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
#'                                    selected = 5e-8,
#'                                    grid     = TRUE),
#'                    sliderTextInput(inputId  = ns("clump_r2"),
#'                                    label    = "r2",
#'                                    choices  = c(0.0001,0.001,0.01,seq(0.1,0.9,0.1),0.9999),
#'                                    selected = 0.001,
#'                                    grid     = TRUE)),
#'             column(6,
#'                    sliderTextInput(inputId  = ns("clump_p2"),
#'                                    label    = "p2",
#'                                    choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
#'                                    selected = 1.0,
#'                                    grid     = TRUE),
#'                    sliderTextInput(inputId  = ns("clump_kb"),
#'                                    label    = "kb",
#'                                    choices  = c(1,seq(0,1000,50)[-1]),
#'                                    selected = 250,
#'                                    grid     = TRUE))
#'           )
#'         )
#'
#'       } else if(input$grouping == "Bayesian") {
#'
#'         controls <- tagList(
#'           fluidRow(
#'             column(6,
#'                    mod_reference_ui(id=ns("reference"))),
#'             column(4,
#'                    actionButton(inputId = ns("finemap"),
#'                                 label   = "Finemap")),
#'             column(1,
#'                    actionButton(inputId = ns("reset"),
#'                                 width   = "40px",
#'                                 label   = "",
#'                                 icon    = icon("rotate-left"))),
#'             tags$style(type='text/css', paste0("#",ns("finemap")," { width:100%; margin-top: 25px;}")),
#'             tags$style(type='text/css', paste0("#",ns("reset")," { width:100%; margin-top: 25px;}")),
#'           ),
#'           fluidRow(
#'             column(3,
#'                    numericInput(inputId  = ns("l_param"),
#'                                 label    = "L",
#'                                 value    = 1,
#'                                 step     = 1)),
#'             column(3,
#'                    numericInput(inputId  = ns("n_param"),
#'                                 label    = "N",
#'                                 value    = NA_real_,
#'                                 step     = 1)),
#'             column(6,
#'                    sliderTextInput(inputId  = ns("finemap_p1"),
#'                                    label    = "p1",
#'                                    choices  = c(5e-8,5e-6,1e-4,0.01,0.05,0.1,0.5,1.0),
#'                                    selected = 1e-4,
#'                                    grid     = TRUE))
#'           ),
#'           fluidRow(
#'             column(3,
#'                    prettyRadioButtons(inputId = ns("source_type"),
#'                                       label   = "Type",
#'                                       choices = c("quant","cc"),
#'                                       selected = "quant")),
#'             column(3,
#'                    numericInput(inputId = ns("sd_y"),
#'                                 label = "sdY",
#'                                 value = NULL,
#'                                 step  = 0.1)),
#'             column(6,
#'                    sliderTextInput(inputId  = ns("coverage"),
#'                                    label    = "Coverage %",
#'                                    choices  = c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0),
#'                                    selected = 0.95,
#'                                    grid     = TRUE))
#'           )
#'         )
#'
#'       }
#'
#'       # return the controls UI elements
#'       return(controls)
#'     })
#'
#'
#'     #==========================================
#'     # reset clumping
#'     #==========================================
#'     session$userData[[ns("reset")]] <- observeEvent(input$reset, {
#'
#'       if(all(c("group","index") %in% names(data_module$data))) {
#'         data_module$data[, c("group","index") := NULL]
#'         tmp <- data.table::copy(data_module$data)
#'         data_module$data <- NULL
#'         data_module$data <- tmp
#'       }
#'
#'     })
#'
#'
#'     #==========================================
#'     # Clump button
#'     #==========================================
#'     session$userData[[ns("clump")]] <- observeEvent(input$clump, {
#'
#'       # check data
#'       if(is.null(data_module$data)) return(NULL)
#'
#'       # reset if previously run (issues with factors being reset)
#'       if(any(c("group","index") %in% names(data_module$data))) {
#'
#'         data_module$data$group <- NULL
#'         data_module$data$index <- NULL
#'
#'       }
#'
#'       # get the reference file
#'       plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
#'                                    chrom    = gene_module$chr,
#'                                    from     = gene_module$start - gene_module$flanks_kb*1000,
#'                                    to       = gene_module$end + gene_module$flanks_kb*1000)
#'
#'
#'       shiny::withProgress(message = 'Clumping data', value = 0, {
#'
#'         shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))
#'
#'         # run clumping
#'         data_module$data <- genepi.utils::clump(gwas      = data_module$data,
#'                                                 p1        = input$clump_p1,
#'                                                 p2        = input$clump_p2,
#'                                                 r2        = input$clump_r2,
#'                                                 kb        = input$clump_kb,
#'                                                 plink2    = get_plink2_exe(),
#'                                                 plink_ref = plink_ref) #|> as.data.frame()
#'
#'         # rename to standard index/group
#'         data.table::setnames(data_module$data, c("index","clump"), c("index","group"))
#'
#'         # warn if no clumps found
#'         if(all(is.na(data_module$data$group))) {
#'           shiny::incProgress(3/4, detail = paste("Failed"))
#'           showNotification("Clumping failed: \nno clumps were found with the provided parameters, or plink failed", type="error")
#'         } else {
#'           shiny::incProgress(3/4, detail = paste("Complete"))
#'         }
#'
#'       })
#'
#'     })
#'
#'
#'     #==========================================
#'     # Finemap button
#'     #==========================================
#'     session$userData[[ns("finemap")]] <- observeEvent(input$finemap, {
#'
#'       # check data
#'       if(is.null(data_module$data)) return(NULL)
#'
#'       # reset if previously run (issues with factors being reset)
#'       if(any(c("group","index") %in% names(data_module$data))) {
#'
#'         data_module$data$group <- NULL
#'         data_module$data$index <- NULL
#'
#'       }
#'
#'       shiny::withProgress(message = 'Finemapping data', value = 0, {
#'
#'         shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))
#'
#'         if(input$l_param > 1) {
#'
#'           # get the reference file
#'           plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
#'                                        chrom    = gene_module$chr,
#'                                        from     = gene_module$start - gene_module$flanks_kb*1000,
#'                                        to       = gene_module$end + gene_module$flanks_kb*1000)
#'
#'           # create the LD matrix for the region variants
#'           dat_ld_obj <- genepi.utils::ld_matrix(dat       = data_module$data,
#'                                                 plink2    = get_plink2_exe(),
#'                                                 plink_ref = plink_ref)
#'
#'           # create coloc dataset
#'           # D1 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
#'           #                          type = input$source_type,
#'           #                          sdY  = input$sd_y,
#'           #                          N    = input$n_param,
#'           #                          L    = input$l_param,
#'           #                          ld   = dat_ld_obj[["ld_mat"]])
#'
#'           results <- susieR::susie_rss(bhat  = dat_ld_obj[["dat"]][["BETA"]],
#'                                        shat  = dat_ld_obj[["dat"]][["SE"]],
#'                                        n     = max(dat_ld_obj[["dat"]][["N"]], na.rm=TRUE),
#'                                        R     = dat_ld_obj[["ld_mat"]],
#'                                        var_y = sdY.est(dat_ld_obj[["dat"]][["SE"]]^2, dat_ld_obj[["dat"]][["EAF"]], max(dat_ld_obj[["dat"]][["N"]], na.rm=TRUE)),
#'                                        L     = input$l_param)
#'                                        # estimate_residual_variance = FALSE,
#'                                        # prior_weights = rep(input$finemap_p1, nrow(dat_ld_obj[["dat"]])),
#'                                        # maf = dat_ld_obj[["dat"]][["EAF"]],
#'                                        # maf_thresh = 0.01)
#'
#' #
#' #           # run SuSie
#' #           results <- coloc::runsusie(D1, coverage=input$coverage, prior_weights=rep(input$finemap_p1, length(D1$snp)))
#'
#'           # update the data with grouping
#'           data_module$data <- calc_credible_set(result      = results,
#'                                                 dat_join_to = data_module$data,
#'                                                 coverage    = input$coverage,
#'                                                 result_type = "coloc::runsusie")
#'
#'           # sensitivity plot
#'           data_module$sensitivity_data <- susieR::kriging_rss(z = dat_ld_obj[["dat"]][["BETA"]]/(dat_ld_obj[["dat"]][["SE"]]),
#'                                                               R = dat_ld_obj[["ld_mat"]],
#'                                                               n = max(dat_ld_obj[["dat"]][["N"]], na.rm=TRUE))$conditional_dist
#'           data_module$sensitivity_data$RSID <- dat_ld_obj[["dat"]][["RSID"]]
#'
#'           data_module$sensitivity_data$outlier <- ifelse(data_module$sensitivity_data$logLR  > 2 &
#'                                                          abs(data_module$sensitivity_data$z) > 2, TRUE, FALSE)
#'
#'         } else {
#'
#'           # create coloc dataset
#'           D1 <- make_coloc_dataset(dat  = data_module$data,
#'                                    type = input$source_type,
#'                                    sdY  = input$sd_y,
#'                                    N    = input$n_param,
#'                                    L    = input$l_param,
#'                                    ld   = NULL)
#'
#'           # if assuming only one significant variant run coloc::finemap.abf()
#'           results <- coloc::finemap.abf(D1, p1=input$finemap_p1)
#'
#'           # update the data with grouping
#'           data_module$data <- calc_credible_set(result      = results,
#'                                                 dat_join_to = data_module$data,
#'                                                 coverage    = input$coverage,
#'                                                 result_type = "coloc::finemap.abf")
#'
#'         }
#'
#'         # warn if no clumps found
#'         if(all(is.na(data_module$data$group))) {
#'           shiny::incProgress(3/4, detail = paste("Failed"))
#'           showNotification("Finemapping failed: \nno groups were found with the provided parameters, or coloc failed", type="error")
#'         } else {
#'           shiny::incProgress(3/4, detail = paste("Complete"))
#'         }
#'
#'       })
#'
#'
#'
#'     })
#'
#'
#'     #==========================================
#'     # Return the the active values for this module
#'     # just indicates whether it is in use/active or
#'     # or not.
#'     #==========================================
#'     return(v)
#'   })
#' }
