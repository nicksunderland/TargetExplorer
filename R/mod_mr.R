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
                               hr(),
                               fluidRow(
                                 column(6, prettyCheckboxGroup(inputId = ns("mr_method"),
                                                               label   = "MR method",
                                                               choices = c("mr_wald_ratio","mr_egger_regression","mr_weighted_median",
                                                                            "mr_ivw","mr_simple_mode", "mr_weighted_mode"),
                                                               selected= c("mr_wald_ratio", "mr_egger_regression", "mr_weighted_median",
                                                                            "mr_ivw", "mr_simple_mode", "mr_weighted_mode"),
                                                               shape   = "curve",
                                                               inline  = FALSE)),
                                 column(6,
                                        prettyCheckboxGroup(inputId = ns("mr_analysis"),
                                                            label   = "Sensitivity",
                                                            choices = c("mr_singlesnp","mr_leaveoneout"),
                                                            selected= c("mr_singlesnp","mr_leaveoneout"),
                                                            shape   = "curve",
                                                            inline  = FALSE),
                                        actionButton(inputId = ns("run_mr"),
                                                     label   = "Run MR"),
                                        hr(),
                                        selectInput(inputId = ns("sens_plot_select"),
                                                    label   = "Plot",
                                                    choices = c("mr_singlesnp","mr_leaveoneout")),
                                        prettyCheckbox(inputId = ns("point_labels"),
                                                       value   = TRUE,
                                                       label   = "Labels",
                                                       shape   = "curve")
                                        )
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
                                                brush    = ns("mr_plot_brush")),
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

    # R CMD checks
    beta.exposure <- beta.outcome <- bp <- RSID <- SNP <- NULL

    #==========================================
    # Data module server for the MR module
    #==========================================
    data_mod     <- mod_data_server(id="data", gene_module=app$modules$gene)
    remove_mod   <- mod_remove_server(id="remove", app=app, parent_id=id, parent_inputs=input)
    exposure_mod <- mod_source_select_server(id="exposure", app=app, source_type=c("GWAS","eQTL","Coloc"), label="Exposure")
    outcome_mod  <- mod_source_select_server(id="outcome",  app=app, source_type=c("GWAS","eQTL","Coloc"), label="Outcome")
    exposure_filter_mod <- mod_source_select_server(id="exposure_filter_by", app=app, source_type=c("GWAS","eQTL","Coloc"), label="filter by:")
    outcome_filter_mod  <- mod_source_select_server(id="outcome_filter_by",  app=app, source_type=c("GWAS","eQTL","Coloc"), label="filter by:")


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
        iv_selection <- c("eqtl","coloc","index")

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

        # exposure data
        shiny::incProgress(1/n, detail = "Formatting exposure")
        exp <- TwoSampleMR::format_data(dat = exposure_mod$data[RSID %in% variants_source_1, ] |> as.data.frame(),
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

        shiny::incProgress(1/n, detail = "Harmonising")
        data_mod$data  <- TwoSampleMR::harmonise_data(exp,out)

        shiny::incProgress(1/n, detail = "MR analysis")
        data_mod$data2 <- TwoSampleMR::mr(data_mod$data)


        # sensitivity if requested
        if("mr_singlesnp" %in% input$mr_analysis) {

          shiny::incProgress(1/n, detail = "Single SNP analysis")
          res_single <- TwoSampleMR::mr_singlesnp(data_mod$data)
          data_mod$plot_single <- TwoSampleMR::mr_forest_plot(res_single)

        }
        if("mr_leaveoneout" %in% input$mr_analysis) {

          shiny::incProgress(1/n, detail = "Leave-one-out SNP analysis")
          res_loo <- TwoSampleMR::mr_leaveoneout(data_mod$data)
          data_mod$plot_loo <- TwoSampleMR::mr_leaveoneout_plot(res_loo)

        }

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

      # create the MR plot
      p   <- TwoSampleMR::mr_scatter_plot(data_mod$data2, data_mod$data)[[1]] +
        theme_classic() +
        theme(legend.position="top")

      # plot labels
      if(input$point_labels) {
        p <- p +
          geom_label_repel(mapping = aes(label = SNP, x = beta.exposure, y = beta.outcome), color="black", show.legend = FALSE)
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

      # clean up the results table
      res <- data_mod$data2
      res$pval <- formatC(res$pval, digits = 3)
      res$id.exposure <- NULL
      res$id.outcome  <- NULL

      return(res)
    })


    #==========================================
    # MR Egger results table
    #==========================================
    output$mr_egg_result <- renderTable({

      # check data
      validate(
        need(!is.null(data_mod$data), "")
      )

      # run analysis
      egg <- TwoSampleMR::mr_pleiotropy_test(data_mod$data)
      if(!is.na(egg$pval)) {
        egg$pval <- formatC(egg$pval, digits = 3)
      }
      egg$id.exposure <- NULL
      egg$id.outcome  <- NULL

      return(egg)
    })



    #==========================================
    # Sensitivity plots
    #==========================================
    output$sens_plot <- renderPlot({

      if(input$sens_plot_select=="mr_singlesnp") {

        # check data
        validate(
          need(!is.null(data_mod$plot_single), paste0('No single SNP analysis data'))
        )
        return(data_mod$plot_single)

      } else if(input$sens_plot_select=="mr_leaveoneout") {

        # check data
        validate(
          need(!is.null(data_mod$plot_loo), paste0('No LOO SNP analysis data'))
        )
        return(data_mod$plot_loo)

      }

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
