#' MR UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_mr_ui <- function(id){
  ns <- NS(id)
  # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
  div(
    id = id,
    sidebarLayout(position = "right",
                  sidebarPanel(width = 3,
                               fluidRow(column(10, p(strong(paste0("Controls [",id,"]")))),
                                        column(2,  mod_remove_ui(id=ns("remove")))),
                               hr(),
                               fluidRow(
                                 column(6, mod_source_select_ui(ns("exposure"))),
                                 column(6, mod_source_select_ui(ns("outcome")))
                               ),
                               fluidRow(
                                 column(6, mod_source_select_ui(ns("exposure_filter_by"))),
                                 column(6, mod_source_select_ui(ns("outcome_filter_by")))
                               ),
                               fluidRow(
                                 shinyjs::disabled(column(6, numericInput(inputId=ns("n_param1"), label="N", value=NA_real_, step=1))),
                                 shinyjs::disabled(column(6, numericInput(inputId=ns("n_param2"), label="N", value=NA_real_, step=1))),
                               ),
                               fluidRow(
                                 column(6,
                                        checkboxGroupInput(inputId  = ns("mr_method"),
                                                           label    = "MR method",
                                                           choices  = c("mr_ivw","mr_egger","mr_weighted_median","mr_weighted_mode","mr_pcgmm"),
                                                           selected = c("mr_ivw","mr_egger"),
                                                           inline   = FALSE)),
                                 column(6,
                                        checkboxInput(inputId=ns("mr_corr"), label="Correlated instruments", value=FALSE),
                                        mod_reference_ui(id=ns("reference")),
                                        shinyjs::disabled(sliderTextInput(inputId=ns("r2thresh"), label="r2 thresh", selected=0.95, choices=c(0.001,0.01,seq(0.1,0.9,by=0.1),0.95,0.98,0.99), grid=TRUE)))
                               ),
                               fluidRow(
                                 column(6, actionButton(inputId=ns("run_mr"), label="Run MR")),
                                 column(6, textInput(inputId=ns("exclude"), label="Exclude", value="", placeholder = "e.g. rs1234; rs4321"))
                               ),
                               fluidRow(
                                 column(6, selectInput(inputId=ns("sens_plot_select"), label="Sensitivity plot", choices=c("Single SNP","Leave-one-out"), selected="Single SNP")),
                                 column(3, checkboxInput(inputId=ns("point_labels"), value=TRUE, label="Labels")),
                                 tags$style(type='text/css', paste0("#",ns("run_mr")," { width:100%; margin-top: 25px;}")),
                                 tags$style(type='text/css', paste0("#",ns("point_labels")," { width:100%; margin-top: 25px;}"))
                               ),
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            fluidRow(
                              column(9, p(strong("Mendelian Randomisation:"))),
                              column(1, p(strong("MR plot:"), style="text-align:right; margin-top: 6px;")),
                              column(1, numericInput(inputId=ns("mr_plot_index"), label=NULL, value=1, step=1, min=1, max=1))
                            ),
                            fluidRow(
                              column(6, plotOutput(outputId = ns("sens_plot"), height="500px")),
                              column(6, plotOutput(outputId = ns("mr_plot"), height="500px", brush=ns("mr_plot_brush")))
                            ),
                            fluidRow(
                              column(12, tableOutput(outputId = ns("mr_result")))
                            ),
                            fluidRow(
                              column(12, tableOutput(outputId = ns("mr_egg_result")))
                            ),
                            fluidRow(
                              column(12, tableOutput(outputId = ns("mr_table")))
                            )
                  ) # main panel end

    ), # sidebar layout end
    hr()
  ) # div end
}


#' MR Server Functions
#' @import ggplot2 ggrepel
#' @importFrom TwoSampleMR mr_scatter_plot mr
#' @noRd
mod_mr_server <- function(id, app){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #==========================================
    # Data module server for the MR module
    #==========================================
    data_mod      <- mod_data_server(id="data", gene_module=app$modules$gene)
    remove_mod    <- mod_remove_server(id="remove", app=app, parent_id=id, parent_inputs=input)
    exposure_mod  <- mod_source_select_server(id="exposure", app=app, source_type=c("GWAS"), label="Exposure", multiple=TRUE)
    outcome_mod   <- mod_source_select_server(id="outcome",  app=app, source_type=c("GWAS"), label="Outcome",  multiple=FALSE)
    reference_mod <- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=FALSE)
    exposure_filter_mod <- mod_source_select_server(id="exposure_filter_by", app=app, source_type=c("GWAS","Coloc"), label="filter by:")
    outcome_filter_mod  <- mod_source_select_server(id="outcome_filter_by",  app=app, source_type=c("GWAS","Coloc"), label="filter by:")


    #==========================================
    # Correlated MR controls
    #==========================================
    session$userData[[ns("mr_corr")]] <- observeEvent(input$mr_corr, {

      shinyjs::toggleState(id="r2thresh",  condition=input$mr_corr)
      shinyjs::toggleState(id="n_param1",  condition=input$mr_corr)
      shinyjs::toggleState(id="n_param2",  condition=input$mr_corr)
      shinyjs::toggleState(id="reference-reference", condition=input$mr_corr)

      current_methods <- input$mr_method
      if(input$mr_corr) {
        updateCheckboxGroupInput(inputId="mr_method",
                                 choices= c("mr_ivw","mr_egger","mr_weighted_median","mr_weighted_mode","mr_pcgmm"),
                                 selected=current_methods)
      } else {
        updateCheckboxGroupInput(inputId="mr_method",
                                 choices= c("mr_ivw","mr_egger","mr_weighted_median","mr_weighted_mode"),
                                 selected=current_methods)
      }

    }, ignoreNULL=TRUE)


    #==========================================
    # Run MR button
    #==========================================
    session$userData[[ns("run_mr")]] <- observeEvent(input$run_mr, {

      # ensure required input data and some methods
      req(exposure_mod$data, outcome_mod$data)
      req(input$mr_method)

      # progress bar
      shiny::withProgress(message = 'Running MR analysis', value = 0, {
        n=8
        shiny::incProgress(1/n, detail = "Starting")

        # possible instrument selection (i.e. is there an `index` column)
        if("index" %in% names(exposure_mod$data)) {
          variants_source_1 <- exposure_mod$data[index==TRUE, RSID]
        } else {
          variants_source_1 <- exposure_mod$data[["RSID"]]
        }

        # apply the filter column (if present) to get the variants
        if("index" %in% names(outcome_mod$data)) {
          variants_source_2 <- outcome_mod$data[index==TRUE, RSID]
        } else {
          variants_source_2 <- outcome_mod$data[["RSID"]]
        }

        # filter the source 1 variants by those in the "filter by/source_1_join" option
        if(!is.null(exposure_filter_mod$data)) {
          if("index" %in% names(exposure_filter_mod$data)) {
            variants_source_1_join <- exposure_filter_mod$data[index==TRUE, RSID]
          } else {
            variants_source_1_join <- exposure_filter_mod$data[["RSID"]]
          }
          # apply to the source 1 variants
          variants_source_1 <- variants_source_1[variants_source_1 %in% variants_source_1_join]
        }

        # filter the source 2 variants by those in the "filter by/source_2_join" option
        if(!is.null(outcome_filter_mod$data)) {
          if("index" %in% names(outcome_filter_mod$data)) {
            variants_source_2_join <- outcome_filter_mod$data[index==TRUE, RSID]
          } else {
            variants_source_2_join <- outcome_filter_mod$data[["RSID"]]
          }
          variants_source_2 <- variants_source_2[variants_source_2 %in% variants_source_2_join]
        }

        # warn if we have filter out all the variants
        if(length(variants_source_1)==0 || length(variants_source_2)==0) {
          showNotification("No variants left when using `filter by:` option(s)", type="warning")
          return(NULL)
        }

        # custom variants to exclude
        exclude_variants <- c("")
        if(!is.null(input$exclude) && input$exclude != "") {
          exclude_variants <- trimws(strsplit(input$exclude,";",fixed=TRUE)[[1]])
          if(!any(exclude_variants %in% exposure_mod$data$RSID)) {
            showNotification("Provided variants were not found in the dataset - recheck input, n.b. must be ';' separated", type="error")
          }
        }

        # if no matching SNPs in outcome then return NULL and warn
        matching_snps = outcome_mod$data$RSID %in% exposure_mod$data$RSID[ exposure_mod$data$RSID %in% variants_source_1 &
                                                                          !exposure_mod$data$RSID %in% exclude_variants    ] &
                        outcome_mod$data$RSID %in% variants_source_2
        if(!any(matching_snps)) {
          showNotification(paste("No matching outcome SNPs for the", nrow(exp), "exposure SNPs provided"), type="error")
          return(NULL)
        }

        # create MR object (harmonise)
        shiny::incProgress(1/n, detail = "Harmonising")
        mr <- MR()
        mr <- mr_load(mr,
                      exposure = exposure_mod$data[RSID %in% variants_source_1 & !RSID %in% exclude_variants, ],
                      outcome  = outcome_mod$data[matching_snps, ],
                      harmonise_strictness = 2)

        # running the external packages may fail
        tryCatch({
          shiny::incProgress(1/n, detail = "MR analysis")

          #-----------------------------------------------
          # Correlated variables, add LD matrix
          if(input$mr_corr) {

            # get the reference file
            if(grepl("1kGv3", reference_mod$ref_path)) {

              plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                           chrom    = app$modules$gene$chr,
                                           from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
                                           to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)
              plink2   <- get_plink2_exe()
              ukbb_ref <- NULL

            } else if(grepl("UKB_LD", reference_mod$ref_path)) {

              # get the UKBB LD file path (downloads if not in cache)
              ukbb_ref_dt <- genepi.utils::download_ukbb_ld(chr           = app$modules$gene$chr,
                                                            bp_start      = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
                                                            bp_end        = app$modules$gene$end + app$modules$gene$flanks_kb*1000,
                                                            ukbb_ld_cache = file.path(reference_mod$ref_path, "cache"))
              ukbb_ref  <- ukbb_ref_dt$root_file
              plink_ref <- NULL
              plink2    <- NULL
            }

            # create the LD matrix for the region variants (flips exposure effect allele to match reference LD reference allele)
            dat_ld_obj <- genepi.utils::ld_matrix(dat = data.frame(RSID = snps_no_alleles(mr),
                                                                   EA   = mr@ea,
                                                                   OA   = mr@oa,
                                                                   EAF  = mr@eafx[,1],
                                                                   BETA = mr@bx[,1]),
                                                  method    = "r",
                                                  plink2    = plink2,
                                                  plink_ref = plink_ref,
                                                  ukbb_ref  = ukbb_ref)

            # set into the object (and harmonises)
            mr <- set_ld_mat(mr,
                             ld_mat = dat_ld_obj$ld_mat,
                             rsid   = dat_ld_obj$dat$RSID,
                             ref    = dat_ld_obj$dat$EA,
                             alt    = dat_ld_obj$dat$OA,
                             align_ref_ea = TRUE,
                             prune_r2_thresh=input$r2thresh)
          }

          #-----------------------------------------------
          # Run MR
          results <- list()

          # inverse (outcome) variance weighted MR
          if("mr_ivw" %in% input$mr_method) { results[['mr_ivw']] <- mr_ivw(mr) }

          # MR Egger
          if("mr_egger" %in% input$mr_method) { results[['mr_egger']] <- mr_egger(mr) }

          # MR weight median
          if("mr_weighted_median" %in% input$mr_method) { results[['mr_weighted_median']] <- mr_weighted_median(mr) }

          # MR weight mode
          if("mr_weighted_mode" %in% input$mr_method) { results[['mr_weighted_mode']] <- mr_weighted_mode(mr) }

          # Principal componenet MR
          if("mr_pcgmm" %in% input$mr_method) { results[['mr_pcgmm']] <- mr_pcgmm(mr) }

          # update the plot index UI
          updateNumericInput(inputId="mr_plot_index", value=1, max=num_exposures(mr))

          # MR results to plotting df
          data_mod$data  <- mr
          data_mod$data2 <- mr_results_to_plotting(results)

        },
        error=function(e) {

          data_mod$data  <- NULL
          data_mod$data2 <- NULL
          showNotification(paste0("Run MR failed: ", e), type="error")
          return(NULL)

        }) # end tryCatch

        shiny::incProgress(1/n, detail = "Complete")
      }) # end progressBar

    })



    #==========================================
    # MR plot
    #==========================================
    output$mr_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(exposure_mod$data), paste0('No exposure data imported')),
        need(!is.null(outcome_mod$data), paste0('No outcome data imported')),
        need(!is.null(data_mod$data) &&
               !is.null(data_mod$data2) &&
               input$mr_plot_index > 0 &&
               input$mr_plot_index <= num_exposures(data_mod$data), "")
      )

      # create plot
      d <-data_mod$data2[id_idx==input$mr_plot_index, ]
      p <- ggplot(data = mr_results_to_plotting(data_mod$data, id_idx=input$mr_plot_index), mapping = aes(x=x, y=y)) +
        geom_errorbar( mapping = aes(ymin=y-yse, ymax=y+yse), width=0, color="grey") +
        geom_errorbarh(mapping = aes(xmin=x-xse, xmax=x+xse), height=0, color="grey") +
        geom_point() +
        geom_abline(data      = d,
                    mapping   = aes(intercept=Intercept, slope=Estimate, color=Method), show.legend=TRUE) +
        ggplot2::scale_colour_manual(values=c("mr_egger"="#1f78b4","mr_weighted_mode"="#b2df8a", "mr_weighted_median"="#33a02c","mr_ivw"="#e31a1c","mr_pcgmm"="#ff7f00", "#fb9a99","#a6cee3", "#fdbf6f","#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
        theme_classic() +
        labs(colour = "MR method",
             x      = paste("SNP effect on", data_mod$data2[id_idx==input$mr_plot_index, Exposure][[1]]),
             y      = paste("SNP effect on", data_mod$data2[id_idx==input$mr_plot_index, Outcome][[1]])) +
        theme(legend.position="top", legend.direction="vertical") +
        guides(colour=ggplot2::guide_legend(ncol=2))

      # # plot labels
      if(input$point_labels) {
        p <- p +
          geom_label_repel(mapping = aes(label=snp))
      }

      return(p)

    })


    #==========================================
    # Point select table - for brushed point
    #==========================================
    output$mr_table <- renderTable({
      tryCatch({
        bp_mr <- brushedPoints(data_mod$data, input$mr_plot_brush)
        if(nrow(bp_mr)>0) return(bp_mr) else return(NULL)
      },
      error = function(e) {
        return(NULL)
      })
    })


    #==========================================
    # MR results table
    #==========================================
    output$mr_result <- renderTable({

      # check data
      validate(
        need(!is.null(data_mod$data2), "")
      )

      d <- data.table::copy(data_mod$data2)[, id_idx:=NULL]

      return(d)
    })



#     #==========================================
#     # Sensitivity plots
#     #==========================================
#     output$sens_plot <- renderPlot({
#
#       # check data
#       validate(
#         need(!is.null(data_mod$data), paste0('No MR data found'))
#       )
#
#       # run and plot single SNP analysis
#       if(input$sens_plot_select=="Single SNP") {
#
#         # run analysis - from MendelianRandmonisation package
#         if(inherits(data_mod$data, "MRInput")) {
#
#           p <- MendelianRandomization::mr_forest(data_mod$data,
#                                                  alpha = 0.05,
#                                                  snp_estimates = TRUE,
#                                                  methods = "ivw",
#                                                  ordered = TRUE)
#
#           # data from TwoSampleMR package
#         } else {
#
#           res_single <- TwoSampleMR::mr_singlesnp(data_mod$data)
#           p <- TwoSampleMR::mr_forest_plot(res_single)
#
#         }
#
#
#         # run and plot LOO SNP analysis
#       } else if(input$sens_plot_select=="Leave-one-out") {
#
#         # run analysis - from MendelianRandmonisation package
#         if(inherits(data_mod$data, "MRInput")) {
#
#           p <- MendelianRandomization::mr_loo(data_mod$data, alpha = 0.05)
#
#           # data from TwoSampleMR package
#         } else {
#
#           res_loo <- TwoSampleMR::mr_leaveoneout(data_mod$data)
#           p <- TwoSampleMR::mr_leaveoneout_plot(res_loo)
#
#         }
#
#       }
#
#       return(p)
#     })




    #==========================================
    # Return the data module to the reactive
    # values in the main app server such that
    # other modules can access the data to use
    # in their own processes.
    #==========================================
    return(data_mod)
  })
}


