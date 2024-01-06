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
                                        checkboxGroupInput(inputId = ns("mr_method"),
                                                           label   = "MR method",
                                                           choices = c("mr_wald_ratio","mr_egger_regression","mr_weighted_median",
                                                                       "mr_ivw","mr_simple_mode", "mr_weighted_mode"),
                                                           selected= c("mr_wald_ratio", "mr_egger_regression", "mr_weighted_median",
                                                                       "mr_ivw", "mr_simple_mode", "mr_weighted_mode"),
                                                           inline  = FALSE)),
                                 column(6,
                                        checkboxGroupInput(inputId = ns("mr_corr_method"),
                                                           label   = "Correlated MR methods",
                                                           choices = c("mr_ivw","mr_egger","mr_pcgmm"),
                                                           selected= NULL,
                                                           inline  = FALSE),
                                        mod_reference_ui(id=ns("reference")))
                               ),
                               fluidRow(
                                 column(6, textInput(inputId=ns("exclude"), label="Exclude", value="", placeholder = "e.g. rs1234; rs4321")),
                                 shinyjs::disabled(column(6, sliderTextInput(inputId=ns("r2thresh"), label="r2 thresh", selected=0.95, choices=c(0.001,0.01,seq(0.1,0.9,by=0.1),0.95,0.98,0.99), grid=TRUE)))
                               ),
                               fluidRow(
                                 column(4, actionButton(inputId=ns("run_mr"), label="Run MR")),
                                 column(5, selectInput(inputId=ns("sens_plot_select"), label="Sensitivity plot", choices=c("Single SNP","Leave-one-out"), selected="Single SNP")),
                                 column(3, checkboxInput(inputId=ns("point_labels"), value=TRUE, label="Labels")),
                                 tags$style(type='text/css', paste0("#",ns("run_mr")," { width:100%; margin-top: 25px;}")),
                                 tags$style(type='text/css', paste0("#",ns("point_labels")," { width:100%; margin-top: 25px;}"))
                               ),
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            p(strong("Mendelian Randomisation:")),
                            fluidRow(
                              column(6,
                                     plotOutput(outputId = ns("sens_plot"),
                                                height   = "500px")
                              ),
                              column(6,
                                     plotOutput(outputId = ns("mr_plot"),
                                                height   = "500px",
                                                brush    = ns("mr_plot_brush")))
                            ),
                            fluidRow(
                              column(12,
                                     tableOutput(outputId = ns("mr_result")),
                                     tableOutput(outputId = ns("mr_egg_result"))
                              )
                            ),
                            tableOutput(outputId = ns("mr_table")),
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
    exposure_mod  <- mod_source_select_server(id="exposure", app=app, source_type=c("GWAS"), label="Exposure")
    outcome_mod   <- mod_source_select_server(id="outcome",  app=app, source_type=c("GWAS"), label="Outcome")
    reference_mod <- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=FALSE)
    exposure_filter_mod <- mod_source_select_server(id="exposure_filter_by", app=app, source_type=c("GWAS","Coloc"), label="filter by:")
    outcome_filter_mod  <- mod_source_select_server(id="outcome_filter_by",  app=app, source_type=c("GWAS","Coloc"), label="filter by:")


    #==========================================
    # Correlated MR controls
    #==========================================
    session$userData[[ns("mr_corr_method")]] <- observeEvent(input$mr_corr_method, {

      shinyjs::toggleState(id="mr_method", condition=is.null(input$mr_corr_method))
      shinyjs::toggleState(id="r2thresh",  condition=!is.null(input$mr_corr_method))
      shinyjs::toggleState(id="n_param1",  condition=!is.null(input$mr_corr_method))
      shinyjs::toggleState(id="n_param2",  condition=!is.null(input$mr_corr_method))
      shinyjs::toggleState(id="reference-reference", condition=!is.null(input$mr_corr_method))

    }, ignoreNULL=FALSE)


    #==========================================
    # Run MR button
    #==========================================
    session$userData[[ns("run_mr")]] <- observeEvent(input$run_mr, {

      # ensure required input data
      req(exposure_mod$data, outcome_mod$data)

      # progress bar
      shiny::withProgress(message = 'Running MR analysis', value = 0, {
        n=8
        shiny::incProgress(1/n, detail = "Starting")

        # possible instrument selection steps
        # in order of how they should be used e.g. if clump and eqtl columns, then use eqtl over clump
        iv_selection <- c("index")

        # get type of data / filter column
        cols_1 <- names(exposure_mod$data)
        cols_2 <- names(outcome_mod$data)
        filter_col_1 <- iv_selection[ which(iv_selection %in% cols_1) ]
        filter_col_2 <- iv_selection[ which(iv_selection %in% cols_2) ]

        # apply the filter column (if present) to get the variants
        if(length(filter_col_1)==0) {
          variants_source_1 <- exposure_mod$data[["RSID"]]
        } else {
          variants_source_1 <- exposure_mod$data[get(filter_col_1[[1]])==TRUE, RSID]
        }

        # apply the filter column (if present) to get the variants
        if(length(filter_col_2)==0) {
          variants_source_2 <- outcome_mod$data[["RSID"]]
        } else {
          variants_source_2 <- outcome_mod$data[get(filter_col_2[[1]])==TRUE, RSID]
        }

        # filter the source 1 variants by those in the "filter by/source_1_join" option
        if(!is.null(exposure_filter_mod$data)) {
          cols_1_join <- names(exposure_filter_mod$data)
          filter_col_1_join <- iv_selection[ which(iv_selection %in% cols_1_join) ]
          if(length(filter_col_1_join)==0) {
            variants_source_1_join <- exposure_filter_mod$data[["RSID"]]
          } else {
            variants_source_1_join <- exposure_filter_mod$data[get(filter_col_1_join[[1]])==TRUE, RSID]
          }
          variants_source_1 <- variants_source_1[variants_source_1 %in% variants_source_1_join]
        }

        # filter the source 2 variants by those in the "filter by/source_2_join" option
        if(!is.null(outcome_filter_mod$data)) {
          cols_2_join <- names(outcome_filter_mod$data)
          filter_col_2_join <- iv_selection[ which(iv_selection %in% cols_2_join) ]
          if(length(filter_col_2_join)==0) {
            variants_source_2_join <- outcome_filter_mod$data[["RSID"]]
          } else {
            variants_source_2_join <- outcome_filter_mod$data[get(filter_col_2_join[[1]])==TRUE, RSID]
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

        # exposure data
        shiny::incProgress(1/n, detail = "Formatting exposure")
        exp <- TwoSampleMR::format_data(dat = exposure_mod$data[RSID %in% variants_source_1 & !RSID %in% exclude_variants, ] |> as.data.frame(),
                                        type = "exposure",
                                        snp_col = "RSID",
                                        beta_col = "BETA",
                                        se_col = "SE",
                                        eaf_col = "EAF",
                                        effect_allele_col = "EA",
                                        other_allele_col = "OA",
                                        pval_col = "P",
                                        chr_col = "CHR",
                                        pos_col = "BP")

        # if no matching SNPs in outcome then return NULL and warn
        matching_snps = outcome_mod$data$RSID %in% exp$SNP &
          outcome_mod$data$RSID %in% variants_source_2
        if(!any(matching_snps)) {
          showNotification(paste("No matching outcome SNPs for the", nrow(exp), "exposure SNPs provided"), type="error")
          return(NULL)
        }

        # format the outcome data
        shiny::incProgress(1/n, detail = "Formatting outcome")
        out <- TwoSampleMR::format_data(dat = outcome_mod$data[matching_snps, ] |> as.data.frame(),
                                        type = "outcome",
                                        snp_col = "RSID",
                                        beta_col = "BETA",
                                        se_col = "SE",
                                        eaf_col = "EAF",
                                        effect_allele_col = "EA",
                                        other_allele_col = "OA",
                                        pval_col = "P",
                                        chr_col = "CHR",
                                        pos_col = "BP")

        # Harmonise
        shiny::incProgress(1/n, detail = "Harmonising")
        data_mod$data  <- TwoSampleMR::harmonise_data(exp,out)
        data_mod$data  <- data_mod$data[data_mod$data$mr_keep==TRUE, ]


        # running the external packages may fail
        tryCatch({

          #-----------------------------------------------
          # Run MR - with independent variables
          shiny::incProgress(1/n, detail = "MR analysis")
          if(is.null(input$mr_corr_method)) {


            data_mod$data2 <- TwoSampleMR::mr(data_mod$data)


          #-----------------------------------------------
          # Run MR - with correlated variables
          } else {

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

            # create the LD matrix for the region variants
            dat_ld_obj <- genepi.utils::ld_matrix(dat       = data_mod$data,
                                                  colmap    = list(RSID = "SNP",
                                                                   EA   = c("effect_allele.exposure","effect_allele.outcome"),
                                                                   OA   = c("other_allele.exposure","other_allele.outcome"),
                                                                   EAF  = c("eaf.exposure","eaf.outcome"),
                                                                   BETA = c("beta.exposure","beta.outcome")),
                                                  method    = "r",
                                                  plink2    = plink2,
                                                  plink_ref = plink_ref,
                                                  ukbb_ref  = ukbb_ref)

            # prune highly correlated variables
            hi_corr_idx       <- prune_hi_corr(dat_ld_obj$ld_mat, thresh=sqrt(input$r2thresh))
            dat_ld_obj$ld_mat <- dat_ld_obj$ld_mat[-hi_corr_idx, -hi_corr_idx]
            dat_ld_obj$dat    <- dat_ld_obj$dat[-hi_corr_idx, ]

            # MRInput object
            data_mod$data <- MendelianRandomization::mr_input(
              bx            = dat_ld_obj$dat$beta.exposure,
              bxse          = dat_ld_obj$dat$se.exposure,
              by            = dat_ld_obj$dat$beta.outcome,
              byse          = dat_ld_obj$dat$se.outcome,
              exposure      = dat_ld_obj$dat$exposure[1],
              outcome       = dat_ld_obj$dat$outcome[1],
              snps          = dat_ld_obj$dat$SNP,
              effect_allele = dat_ld_obj$dat$effect_allele.exposure,
              other_allele  = dat_ld_obj$dat$other_allele.exposure,
              eaf           = dat_ld_obj$dat$eaf.exposure,
              correlation   = dat_ld_obj$ld_mat
            )

            # Run MR methods
            results <- list()
            if("mr_ivw" %in% input$mr_corr_method) {

              method <- "mr_ivw"

              res <- do.call(getFromNamespace(method, "MendelianRandomization"), list(object=data_mod$data, correl=TRUE))
              r <- data.table::data.table(nSNPs     = res@SNPs,
                                          `High R2` = length(hi_corr_idx),
                                          Estimate  = res@Estimate,
                                          SE        = res@StdError,
                                          pval      = res@Pvalue,
                                          Fstat     = res@Fstat,
                                          Intercept = 0)
              results[[method]] <- r

            }

            if("mr_egger" %in% input$mr_corr_method) {

              method <- "mr_egger"

              res <- do.call(getFromNamespace(method, "MendelianRandomization"), list(object=data_mod$data, correl=TRUE))
              r <- data.table::data.table(nSNPs     = res@SNPs,
                                          `High R2` = length(hi_corr_idx),
                                          Estimate  = res@Estimate,
                                          SE        = res@StdError.Est,
                                          pval      = res@Pvalue.Est,
                                          Intercept = res@Intercept,
                                          `SE-int`= res@StdError.Int,
                                          `P-value-int`=res@Pvalue.Int)
              results[[method]] <- r

            }

            if("mr_pcgmm" %in% input$mr_corr_method) {

              method <- "mr_pcgmm"

              # require sample size
              if(is.na(input$n_param1) && !"N" %in% names(exposure_mod$data)) {
                stop("Exposure sample size `N` must be provided for `mr_pcgmm`")
              } else if(is.na(input$n_param1)) {
                nx <- max(exposure_mod$data$N, na.rm=TRUE)
                updateNumericInput(inputId="n_param1", value = nx)
              } else {
                nx <- input$n_param1
              }
              if(is.na(input$n_param2) && !"N" %in% names(outcome_mod$data)) {
                stop("Outcome sample size `N` must be provided for `mr_pcgmm`")
              } else if(is.na(input$n_param2)) {
                ny <- max(outcome_mod$data$N, na.rm=TRUE)
                updateNumericInput(inputId="n_param2", value = ny)
              } else {
                ny <- input$n_param2
              }

              res <- do.call(getFromNamespace(method, "MendelianRandomization"), list(object=data_mod$data, nx=nx, ny=ny))
              r <- data.table::data.table(nSNPs     = nrow(res@Correlation),
                                          `High R2` = length(hi_corr_idx),
                                          Estimate  = res@Estimate,
                                          SE        = res@StdError,
                                          pval      = res@Pvalue,
                                          Intercept = 0,
                                          Fstat     = res@Fstat,
                                          Overdispersion=res@Overdispersion,
                                          PCs       = res@PCs)
              results[[method]] <- r

            }

            # set results
            results <- data.table::rbindlist(results, idcol="Method", fill=TRUE)
            class(results) <- append("MendelianRandomization", class(results))
            data_mod$data2 <- results

          }

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
        need(!is.null(data_mod$data) && !is.null(data_mod$data2), paste0('No MR data found, have you clicked `Run MR`?'))
      )

      # create the MR plot - from MendelianRandmonisation package
      if(inherits(data_mod$data, "MRInput")) {

        plot_data <- data.table::data.table(SNP           = paste0(data_mod$data@snps,"_",data_mod$data@other_allele,"_",data_mod$data@effect_allele),
                                            exposure.beta = data_mod$data@betaX,
                                            outcome.beta  = data_mod$data@betaY,
                                            outcome.se    = data_mod$data@betaXse,
                                            exposure.se   = data_mod$data@betaYse)

        p <- ggplot(data    = plot_data,
                    mapping = aes(x = exposure.beta, y = outcome.beta)) +
          geom_errorbar(mapping = aes(ymin=outcome.beta-(1.96*outcome.se), ymax=outcome.beta+(1.96*outcome.se)), width=0, color="grey") +
          geom_errorbar(mapping = aes(xmin=exposure.beta-(1.96*exposure.se), xmax=exposure.beta+(1.96*exposure.se)), width=0, color="grey") +
          geom_point() +
          geom_abline(data      = data_mod$data2,
                      mapping   = aes(intercept=Intercept, slope=Estimate, color=Method)) +
          theme_classic() +
          theme(legend.position="top")

        # plot labels
        if(input$point_labels) {
          p <- p +
            geom_label_repel(mapping = aes(label=SNP))
        }

        # data from TwoSampleMR package
      } else {

        p   <- TwoSampleMR::mr_scatter_plot(data_mod$data2, data_mod$data)[[1]] +
          theme_classic() +
          theme(legend.position="top")

        # plot labels
        if(input$point_labels) {
          p <- p +
            geom_label_repel(mapping = aes(label = SNP, x = beta.exposure, y = beta.outcome), color="black", show.legend = FALSE)
        }

      }


      # return plot
      return(p)
    })


    #==========================================
    # Point select table - for brushed point
    #==========================================
    output$mr_table <- renderTable({
      tryCatch({
        sign_adj_data <- data_mod$data
        sign_adj_data$beta.exposure <- abs(sign_adj_data$beta.exposure)
        bp_mr <- brushedPoints(sign_adj_data, input$mr_plot_brush)
        if(nrow(bp)>0) return(bp) else return(NULL)
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

      res <- data.table::as.data.table(data_mod$data2)
      res[, pval := formatC(pval, digits = 3)]

      return(res)
    })


    #==========================================
    # MR Egger results table
    #==========================================
    output$mr_egg_result <- renderTable({

      # check data
      validate(
        need(!is.null(data_mod$data2), "")
      )

      # run analysis - from MendelianRandmonisation package
      if(inherits(data_mod$data2, "MendelianRandomization")) {

        egg <- NULL # result lives in the main table

        # data from TwoSampleMR package
      } else {

        egg <- TwoSampleMR::mr_pleiotropy_test(data_mod$data)
        if(!is.na(egg$pval)) {
          egg$pval <- formatC(egg$pval, digits = 3)
        }
        egg$id.exposure <- NULL
        egg$id.outcome  <- NULL

      }

      return(egg)
    })



    #==========================================
    # Sensitivity plots
    #==========================================
    output$sens_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(data_mod$data), paste0('No MR data found'))
      )

      # run and plot single SNP analysis
      if(input$sens_plot_select=="Single SNP") {

        # run analysis - from MendelianRandmonisation package
        if(inherits(data_mod$data, "MRInput")) {

          p <- MendelianRandomization::mr_forest(data_mod$data,
                                                 alpha = 0.05,
                                                 snp_estimates = TRUE,
                                                 methods = "ivw",
                                                 ordered = TRUE)

          # data from TwoSampleMR package
        } else {

          res_single <- TwoSampleMR::mr_singlesnp(data_mod$data)
          p <- TwoSampleMR::mr_forest_plot(res_single)

        }


        # run and plot LOO SNP analysis
      } else if(input$sens_plot_select=="Leave-one-out") {

        # run analysis - from MendelianRandmonisation package
        if(inherits(data_mod$data, "MRInput")) {

          p <- MendelianRandomization::mr_loo(data_mod$data, alpha = 0.05)

          # data from TwoSampleMR package
        } else {

          res_loo <- TwoSampleMR::mr_leaveoneout(data_mod$data)
          p <- TwoSampleMR::mr_leaveoneout_plot(res_loo)

        }

      }

      return(p)
    })




    #==========================================
    # Return the data module to the reactive
    # values in the main app server such that
    # other modules can access the data to use
    # in their own processes.
    #==========================================
    return(data_mod)
  })
}



prune_hi_corr <- function(ld_mat, thresh=sqrt(0.95), seed=2024) {

  # https://wellcomeopenresearch.org/articles/8-449
  # threshold set to r2>0.95
  set.seed(seed)        # for reproducibility
  omit = NULL           # set up list of variants to be omitted
  rho.upper = ld_mat    # correlation matrix
  rho.upper[lower.tri(ld_mat, diag=TRUE)] <- 0
  # only consider upper triangle of correlations
  j=1                   # set counter to 1

  while (max(abs(rho.upper), na.rm=TRUE) > thresh) {
    omit[j] = ifelse(rbinom(1, 1, 0.5)==1,
                     which.max(apply(abs(rho.upper), 1, max, na.rm=TRUE)),
                     which.max(apply(abs(rho.upper), 2, max, na.rm=TRUE)))
    # find the highest correlation value
    # select either the row or column at random
    # add this to the list of omitted variants
    rho.upper[omit[j],] <- 0 # set the correlations in this row to zero
    rho.upper[,omit[j]] <- 0 # set the correlations in this column to zero
    # (to avoid selecting the same variant again)
    j=j+1                # increment the counter
  }  # stop when no more pairwise correlations exceed threshold

  return(omit)
}

