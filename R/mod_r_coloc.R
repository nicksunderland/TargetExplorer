#' r_coloc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#' @importFrom shinyalert shinyalert
#'
mod_r_coloc_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6,
             uiOutput(outputId = ns("function_select"))),
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
    uiOutput(outputId = ns("function_controls"))
  ) # tagList end
}


#' r_coloc Server Functions
#' @noRd
mod_r_coloc_server <- function(id,
                               gene_module,
                               source_1_module,
                               source_2_module = NULL,
                               data_module     = NULL,
                               parent_ui       = NULL,
                               functions = c("finemap.abf","finemap.signals","coloc.abf","coloc.signals","coloc.susie")) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #==========================================
    # Reactive values
    #==========================================
    run_flag <- reactiveVal(FALSE)


    #==========================================
    # Reference (sub)module server for the r-coloc module
    #==========================================
    reference_mod <- NULL


    #==========================================
    # Run button
    #==========================================
    session$userData[[ns("run")]] <- observeEvent(input$run, {

      # check data
      req(source_1_module$data)

      # analysing for multiple signals can take ages with lots of points, confirm run first
      if(input$function_select %in% c("finemap.signals","coloc.signals","coloc.susie") && nrow(source_1_module$data) > 2000) {

        shinyalert(
          inputId = "shinyalert",
          title = "Processing time warning",
          text = paste0(nrow(source_1_module$data), " data points currently loaded, this may take some time. \n",
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

      # functions might fail
      tryCatch({

        shiny::withProgress(message = 'r-coloc package', value = 0, {

          #--------------------------------------------
          # finemap.abf
          if(input$function_select == "finemap.abf")
          {
            shiny::incProgress(1/4, detail = "finemap.abf")

            # reset if previously run
            if(any(c("group","index") %in% names(source_1_module$data))) {

              source_1_module$data$group <- NULL
              source_1_module$data$index <- NULL

            }

            # create coloc dataset
            D1 <- make_coloc_dataset(dat  = source_1_module$data,
                                     type = input$source_type1,
                                     sdY  = input$sd_y1,
                                     N    = input$n_param1)

            # if assuming only one significant variant run coloc::finemap.abf()
            results <- coloc::finemap.abf(dataset = D1, p1 = input$p1)

            # store top 5 as a summary table
            data_module[["summary_table"]] <- results[order(-results$SNP.PP), ] |> head(5)

            # update the data with grouping
            updated_dat <- source_1_module$data[RSID == results$snp[which.max(results$SNP.PP)], c("index","group") := list(TRUE, factor(1))]
            source_1_module$data <- NULL # force memory location change for reactives
            source_1_module$data <- updated_dat

          #--------------------------------------------
          # finemap.signals
          }
          else if(input$function_select == "finemap.signals")
          {
            shiny::incProgress(1/4, detail = "finemap.signals")

            # requirements and checks
            req(source_1_module$data, reference_mod$ref_path)

            # reset if previously run
            if(any(c("group","index") %in% names(source_1_module$data))) {

              source_1_module$data$group <- NULL
              source_1_module$data$index <- NULL

            }

            # get LD matrix and harmonise data
            dat_ld_obj <- get_harmonised_ld(dat      = source_1_module$data,
                                            ref_path = reference_mod$ref_path,
                                            chr      = gene_module$chr,
                                            from     = gene_module$start - gene_module$flanks_kb*1000,
                                            to       = gene_module$end + gene_module$flanks_kb*1000,
                                            method   = "r")

            # store the r corr matrix
            data_module$ld_matrix_r <- dat_ld_obj[["ld_mat"]]

            # reset the r2 matrix
            data_module$ld_matrix_r2 <- NULL

            # make and check the dataset
            D1 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
                                     type = input$source_type1,
                                     sdY  = input$sd_y1,
                                     S    = input$s_param1,
                                     N    = input$n_param1,
                                     ld   = dat_ld_obj[["ld_mat"]])

            # run finemap.signals
            results <- coloc::finemap.signals(D         = D1,
                                              LD        = D1$LD,
                                              method    = input$method,
                                              r2thr     = 0.01,
                                              sigsnps   = NULL,
                                              pthr      = input$p1,
                                              maxhits   = ifelse(input$method=="single",1,input$l_param),
                                              return.pp = FALSE)

            # store top 5 as a summary table
            data_module[["summary_table"]] <- data.frame(t(results))

            # update the data with grouping
            updated_dat <- source_1_module$data[RSID %in% dat_ld_obj$dat[LD_RSID_allele %in% names(results), RSID], index := TRUE]
            updated_dat[index==TRUE, group := factor(1:.N)]
            source_1_module$data <- NULL # force memory location change for reactives
            source_1_module$data <- updated_dat

            # store data for the kriging_rss plot
            if(!is.null(parent_ui) && input$method=="cond" && parent_ui$sensitivity_plot()=="Kriging plot") {
              source_1_module$kriging_rss         <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), D1$LD, n=D1$N)$conditional_dist
              source_1_module$kriging_rss$RSID    <- D1$snp
              source_1_module$kriging_rss$outlier <- ifelse(source_1_module$kriging_rss$logLR  > 2 &
                                                            abs(source_1_module$kriging_rss$z) > 2, TRUE, FALSE)
            } else {
              source_1_module$kriging_rss <- NULL
            }


          #--------------------------------------------
          # coloc.abf
          }
          else if(input$function_select == "coloc.abf")
          {
            shiny::incProgress(1/4, detail = "coloc.abf")

            # require data from 2 input datasets for colocalisation
            req(source_1_module$data, source_2_module$data)

            # reset if previously run
            if(!is.null(data_module$data)) {

              data_module$data <- NULL

            }

            # harmonise the datasets
            harm <- genepi.utils::harmonise(gwas1 = source_1_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P,N)],
                                            gwas2 = source_2_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P,N)],
                                            gwas1_trait = "d1",
                                            gwas2_trait = "d2",
                                            merge = c("SNP"="SNP"))
            harm <- harm[keep==TRUE, ]

            # make and check the dataset 1
            D1 <- make_coloc_dataset(dat  = harm[, list(RSID=RSID_d1,BP=BP_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
                                     type = input$source_1_type,
                                     sdY  = input$sd_y1,
                                     N    = input$n_param1,
                                     S    = input$s_param1)

            # make and check the dataset 2
            D2 <- make_coloc_dataset(dat  = harm[, list(RSID=RSID_d2,BP=BP_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
                                     type = input$source_2_type,
                                     sdY  = input$sd_y2,
                                     N    = input$n_param2,
                                     S    = input$s_param2)

            # run the coloc abf function
            res <- coloc::coloc.abf(dataset1 = D1,
                                    dataset2 = D2,
                                    p1       = input$p1,
                                    p2       = input$p2,
                                    p12      = input$p12)

            # store the summary table
            data_module[["summary_table"]] <- t(as.data.frame(res$summary))

            # create the coloc data.table
            coloc_table <- coloc_to_plotting(res, h4=input$h4, coverage=input$coverage)
            coloc_table[harm, c("BP_1","nlog10P_1") := list(BP_d1,nlog10P_d1), on=c("RSID_1"="RSID_d1")]
            coloc_table[harm, c("BP_2","nlog10P_2") := list(BP_d2,nlog10P_d2), on=c("RSID_2"="RSID_d2")]

            # set
            coloc_table[, RSID := RSID_1] # need to maintain an RSID column, e.g. for MR module to filter on
            data_module$data <- coloc_table

            # store the sensitivity plots - might fail with unusual prior probability combinations; so do last.
            data_module[["coloc_prob_prior"]][[1]] <- genepi.utils::plot_coloc_probabilities(res, rule=paste0("H4 > ", input$h4), type="prior")
            data_module[["coloc_prob_post"]][[1]]  <- genepi.utils::plot_coloc_probabilities(res, rule=paste0("H4 > ", input$h4), type="posterior")


          #--------------------------------------------
          # coloc.signals
          }
          else if(input$function_select == "coloc.signals")
          {
            shiny::incProgress(1/4, detail = "coloc.signals")

            # harmonise the datasets
            harm <- genepi.utils::harmonise(gwas1 = source_1_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P,N)],
                                            gwas2 = source_2_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P,N)],
                                            gwas1_trait = "d1",
                                            gwas2_trait = "d2",
                                            merge = c("SNP"="SNP"))
            harm <- harm[keep==TRUE, ]

            # colmap for getting harmonised LD
            colmap=list(RSID = c("RSID_d1", "RSID_d1"),
                        EA   = c("EA_d1","EA_d2"),
                        OA   = c("OA_d1","OA_d2"),
                        BETA = c("BETA_d1","BETA_d2"),
                        EAF  = c("EAF_d1","EAF_d2"))

            # get LD matrix and harmonise data
            dat_ld_obj <- get_harmonised_ld(dat      = harm,
                                            colmap   = colmap,
                                            ref_path = reference_mod$ref_path,
                                            chr      = gene_module$chr,
                                            from     = gene_module$start - gene_module$flanks_kb*1000,
                                            to       = gene_module$end + gene_module$flanks_kb*1000,
                                            method   = "r")

            # make and check the dataset 1
            D1 <- make_coloc_dataset(dat  = dat_ld_obj[['dat']][, list(LD_RSID_allele,RSID=RSID_d1,BP=BP_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
                                     type = input$source_1_type,
                                     sdY  = input$sd_y1,
                                     N    = input$n_param1,
                                     S    = input$s_param1,
                                     ld   = dat_ld_obj[['ld_mat']])

            # make and check the dataset 2
            D2 <- make_coloc_dataset(dat  = dat_ld_obj[['dat']][, list(LD_RSID_allele,RSID=RSID_d2,BP=BP_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
                                     type = input$source_2_type,
                                     sdY  = input$sd_y2,
                                     N    = input$n_param2,
                                     S    = input$s_param2,
                                     ld   = dat_ld_obj[['ld_mat']])



            # run the coloc signals function
            res <- coloc::coloc.signals(dataset1 = D1,
                                        dataset2 = D2,
                                        method   = input$method, #c("single", "cond", "mask"),
                                        mode     = input$mode, #c("iterative", "allbutone"),
                                        p1       = input$p1, #1e-04,
                                        p2       = input$p2, #1e-04,
                                        p12      = input$p12, #1e-5
                                        maxhits  = input$l_param, # 3
                                        r2thr    = input$r2thr, #0.01,
                                        pthr     = input$p_secondary) # 1e-06 #)

            # store the summary table
            data_module[["summary_table"]] <- as.data.frame(res$summary)

            # store the sensitivity plots
            for(row in 1:nrow(res$summary)) {

              data_module[["coloc_prob_prior"]][[row]] <- genepi.utils::plot_coloc_probabilities(res, rule=paste0("H4 > ", input$h4), type="prior", row=row)
              data_module[["coloc_prob_post"]][[row]]  <- genepi.utils::plot_coloc_probabilities(res, rule=paste0("H4 > ", input$h4), type="posterior", row=row)

            }

            # create the coloc data.table
            coloc_table <- coloc_to_plotting(res, h4=input$h4, coverage=input$coverage)
            coloc_table[dat_ld_obj[['dat']], c("BP_1","nlog10P_1") := list(BP_d1,nlog10P_d1), on=c("RSID_1"="LD_RSID_allele")]
            coloc_table[dat_ld_obj[['dat']], c("BP_2","nlog10P_2") := list(BP_d2,nlog10P_d2), on=c("RSID_2"="LD_RSID_allele")]

            # set
            coloc_table[, RSID := RSID_1] # need to maintain an RSID column, e.g. for MR module to filter on
            data_module$data <- coloc_table

            # store data for the kriging_rss plots
            if(!is.null(parent_ui) && input$method=="cond" && parent_ui$sensitivity_plot()=="Kriging plot") {

              data_module[["kriging_rss"]]         <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), D1$LD, n=D1$N)$conditional_dist
              data_module[["kriging_rss"]]$RSID    <- D1$snp
              data_module[["kriging_rss"]]$outlier <- ifelse(data_module$kriging_rss$logLR  > 2 &
                                                              abs(data_module$kriging_rss$z) > 2, TRUE, FALSE)
              data_module[["kriging_rss2"]]         <- susieR::kriging_rss(D2$beta/(D2$varbeta ^ 0.5), D1$LD, n=D2$N)$conditional_dist
              data_module[["kriging_rss2"]]$RSID    <- D2$snp
              data_module[["kriging_rss2"]]$outlier <- ifelse(data_module$kriging_rss2$logLR  > 2 &
                                                               abs(data_module$kriging_rss2$z) > 2, TRUE, FALSE)

            } else {

              data_module[["kriging_rss"]]  <- NULL
              data_module[["kriging_rss2"]] <- NULL

            }

          #--------------------------------------------
          # coloc.susie
          }
          else if(input$function_select == "coloc.susie")
          {
            shiny::incProgress(1/4, detail = "coloc.susie")

            # harmonise the datasets
            harm <- genepi.utils::harmonise(gwas1 = source_1_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P,N)],
                                            gwas2 = source_2_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P,N)],
                                            gwas1_trait = "d1",
                                            gwas2_trait = "d2",
                                            merge = c("SNP"="SNP"))
            harm <- harm[keep==TRUE, ]

            # colmap for getting harmonised LD
            colmap=list(RSID = c("RSID_d1", "RSID_d1"),
                        EA   = c("EA_d1","EA_d2"),
                        OA   = c("OA_d1","OA_d2"),
                        BETA = c("BETA_d1","BETA_d2"),
                        EAF  = c("EAF_d1","EAF_d2"))

            # get LD matrix and harmonise data
            dat_ld_obj <- get_harmonised_ld(dat      = harm,
                                            colmap   = colmap,
                                            ref_path = reference_mod$ref_path,
                                            chr      = gene_module$chr,
                                            from     = gene_module$start - gene_module$flanks_kb*1000,
                                            to       = gene_module$end + gene_module$flanks_kb*1000,
                                            method   = "r")

            # make and check the dataset 1
            D1 <- make_coloc_dataset(dat  = dat_ld_obj[['dat']][, list(LD_RSID_allele,RSID=RSID_d1,BP=BP_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
                                     type = input$source_1_type,
                                     sdY  = input$sd_y1,
                                     N    = input$n_param1,
                                     S    = input$s_param1,
                                     ld   = dat_ld_obj[['ld_mat']])

            # make and check the dataset 2
            D2 <- make_coloc_dataset(dat  = dat_ld_obj[['dat']][, list(LD_RSID_allele,RSID=RSID_d2,BP=BP_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
                                     type = input$source_2_type,
                                     sdY  = input$sd_y2,
                                     N    = input$n_param2,
                                     S    = input$s_param2,
                                     ld   = dat_ld_obj[['ld_mat']])

            # run the coloc signals function
            res <- coloc::coloc.susie(dataset1   = D1,
                                      dataset2   = D2,
                                      p1         = input$p1, #1e-04,
                                      p2         = input$p2, #1e-04,
                                      p12        = input$p12, #1e-5
                                      susie.args = list(L          = input$l_param,
                                                        coverage   = input$coverage,
                                                        maf_thresh = input$maf_thresh))

            # store the summary table
            data_module[["summary_table"]] <- as.data.frame(res$summary)

            # store the sensitivity plots
            for(row in 1:nrow(res$summary)) {

              data_module[["coloc_prob_prior"]][[row]] <- genepi.utils::plot_coloc_probabilities(res, rule=paste0("H4 > ", input$h4), type="prior", row=row)
              data_module[["coloc_prob_post"]][[row]]  <- genepi.utils::plot_coloc_probabilities(res, rule=paste0("H4 > ", input$h4), type="posterior", row=row)

            }

            # create the coloc data.table
            coloc_table <- coloc_to_plotting(res, h4=input$h4, coverage=input$coverage)
            coloc_table[dat_ld_obj[['dat']], c("BP_1","nlog10P_1") := list(BP_d1,nlog10P_d1), on=c("RSID_1"="LD_RSID_allele")]
            coloc_table[dat_ld_obj[['dat']], c("BP_2","nlog10P_2") := list(BP_d2,nlog10P_d2), on=c("RSID_2"="LD_RSID_allele")]

            # set
            coloc_table[, RSID := RSID_1] # need to maintain an RSID column, e.g. for MR module to filter on
            data_module$data <- coloc_table

            # store data for the kriging_rss plots
            if(!is.null(parent_ui) && parent_ui$sensitivity_plot()=="Kriging plot") {

              data_module[["kriging_rss"]]         <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), D1$LD, n=D1$N)$conditional_dist
              data_module[["kriging_rss"]]$RSID    <- D1$snp
              data_module[["kriging_rss"]]$outlier <- ifelse(data_module$kriging_rss$logLR  > 2 &
                                                               abs(data_module$kriging_rss$z) > 2, TRUE, FALSE)
              data_module[["kriging_rss2"]]         <- susieR::kriging_rss(D2$beta/(D2$varbeta ^ 0.5), D1$LD, n=D2$N)$conditional_dist
              data_module[["kriging_rss2"]]$RSID    <- D2$snp
              data_module[["kriging_rss2"]]$outlier <- ifelse(data_module$kriging_rss2$logLR  > 2 &
                                                                abs(data_module$kriging_rss2$z) > 2, TRUE, FALSE)

            }

          } # end different function type runs


          # warn if no groups found
          if(all(is.na(c(data_module$data$group, source_1_module$data$group)))) {
            shiny::incProgress(3/4, detail = paste("Failed"))
            showNotification("Finemapping failed: \nno groups were found with the provided parameters, or coloc failed", type="error")
          } else {
            shiny::incProgress(3/4, detail = paste("Complete"))
          }

        }) # end withProgress

      },
      error=function(e) {

        showNotification(paste0(input$function_select, " failed - ", e), type="error", duration = 30)
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
    # Function select - which functions to allow
    #==========================================
    output$function_select <- renderUI({

      selectInput(inputId=ns("function_select"), label="Function", choices=functions, selected=functions[[1]])

    })


    #==========================================
    # Function UI - which functionUI do we need
    #==========================================
    output$function_controls <- renderUI({

      req(input$function_select)


      #--------------------------------------------
      # finemap.abf
      if(input$function_select == "finemap.abf")
      {

        controls <- tagList(
          fluidRow(
            column(3, prettyRadioButtons(inputId=ns("source_type1"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("s_param1"), label="S",   value=NA_real_, step=0.1)),
            column(3, numericInput(inputId=ns("sd_y1"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, numericInput(inputId=ns("n_param1"), label="N",   value=NA_real_, step=1)),
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p1"), label="p1", selected=1e-4, choices=c(5e-8,5e-6,1e-4,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("coverage"), label="Coverage %", selected=0.95, choices=c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0), grid=TRUE))
          )
        ) # end tagList


      #--------------------------------------------
      # finemap.signals
      }
      else if(input$function_select == "finemap.signals")
      {

        controls <- tagList(
          fluidRow(
            column(6, mod_reference_ui(id=ns("reference"))),
            column(3, prettyRadioButtons(inputId=ns("source_type1"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, prettyRadioButtons(inputId=ns("method"), label="Method", choices=c("mask","cond","single"), selected="mask", inline=FALSE))
          ),
          fluidRow(
            column(3, numericInput(inputId=ns("l_param"),  label="L",   value=5,        step=1)),
            column(3, numericInput(inputId=ns("s_param1"), label="S",   value=NA_real_, step=0.1)),
            column(3, numericInput(inputId=ns("sd_y1"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, numericInput(inputId=ns("n_param1"), label="N",   value=NA_real_, step=1)),
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p1"), label="p1", selected=1e-4, choices=c(5e-8,5e-6,1e-4,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("coverage"), label="Coverage %", selected=0.95, choices=c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0), grid=TRUE))
          )
        ) # end tagList

        reference_mod <<- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=FALSE)


      #--------------------------------------------
      # coloc.abf
      }
      else if(input$function_select == "coloc.abf")
      {

        controls <- tagList(
          # Type & N-sample size
          fluidRow(
            #------- 1
            column(3, prettyRadioButtons(inputId=ns("source_type1"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("n_param1"), label="N",   value=NA_real_, step=1)),
            #------- 2
            column(3, prettyRadioButtons(inputId=ns("source_type2"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("n_param2"), label="N",   value=NA_real_, step=1))
          ),
          # SdY and S-proportion of sample that are cases
          fluidRow(
            #------- 1
            column(3, numericInput(inputId=ns("sd_y1"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, shinyjs::disabled(numericInput(inputId=ns("s_param1"), label="S",   value=NA_real_, step=0.1))),
            #------- 2
            column(3, numericInput(inputId=ns("sd_y2"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, shinyjs::disabled(numericInput(inputId=ns("s_param2"), label="S",   value=NA_real_, step=0.1))),
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p1"), label="p1", selected=1e-4, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("p2"), label="p2", selected=1e-4, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p12"), label="p12", selected=1e-5, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("h4"), label="H4 decision", selected=0.5, choices=seq(0,1,by=0.1), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("coverage"), label="Coverage %", selected=0.95, choices=c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0), grid=TRUE))
          )

        ) # end tagList


      #--------------------------------------------
      # coloc.signals
      }
      else if(input$function_select == "coloc.signals")
      {

        controls <- tagList(
          # Type & N-sample size
          fluidRow(
            #------- 1
            column(3, prettyRadioButtons(inputId=ns("source_type1"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("n_param1"), label="N",   value=NA_real_, step=1)),
            #------- 2
            column(3, prettyRadioButtons(inputId=ns("source_type2"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("n_param2"), label="N",   value=NA_real_, step=1))
          ),
          # SdY and S-proportion of sample that are cases
          fluidRow(
            #------- 1
            column(3, numericInput(inputId=ns("sd_y1"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, shinyjs::disabled(numericInput(inputId=ns("s_param1"), label="S",   value=NA_real_, step=0.1))),
            #------- 2
            column(3, numericInput(inputId=ns("sd_y2"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, shinyjs::disabled(numericInput(inputId=ns("s_param2"), label="S",   value=NA_real_, step=0.1))),
          ),
          fluidRow(
            column(6, mod_reference_ui(id=ns("reference"))),
            column(3, shinyjs::disabled(prettyRadioButtons(inputId=ns("mode"), label="Mode", choices=c("iterative","allbutone"), selected="iterative", inline=FALSE))),
            column(3, prettyRadioButtons(inputId=ns("method"), label="Method", choices=c("mask","cond","single"), selected="mask", inline=FALSE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p1"), label="p1", selected=1e-4, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("p2"), label="p2", selected=1e-4, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p12"), label="p12", selected=1e-5, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("h4"), label="H4 decision", selected=0.5, choices=seq(0,1,by=0.1), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("r2thr"), label="R2 thresh", selected=0.01, choices=c(0.0001,0.001,0.01,seq(0.1,0.9,0.1),0.9999), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("p_secondary"), label="p_secondary", selected=1e-6, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE))
          ),
          fluidRow(
            column(6, numericInput(inputId=ns("l_param"), label="L",   value=3,        step=1)),
            column(6, sliderTextInput(inputId=ns("coverage"), label="Coverage %", selected=0.95, choices=c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0), grid=TRUE))
          )
        ) # end tagList

        reference_mod <<- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=FALSE)


      #--------------------------------------------
      # coloc.susie
      }
      else if(input$function_select == "coloc.susie")
      {

        controls <- tagList(
          # Type & N-sample size
          fluidRow(
            #------- 1
            column(3, prettyRadioButtons(inputId=ns("source_type1"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("n_param1"), label="N",   value=NA_real_, step=1)),
            #------- 2
            column(3, prettyRadioButtons(inputId=ns("source_type2"), label="Type", choices=c("quant","cc"), selected="quant", inline=FALSE)),
            column(3, numericInput(inputId=ns("n_param2"), label="N",   value=NA_real_, step=1))
          ),
          # SdY and S-proportion of sample that are cases
          fluidRow(
            #------- 1
            column(3, numericInput(inputId=ns("sd_y1"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, shinyjs::disabled(numericInput(inputId=ns("s_param1"), label="S",   value=NA_real_, step=0.1))),
            #------- 2
            column(3, numericInput(inputId=ns("sd_y2"),    label="sdY", value=NA_real_, step=0.1)),
            column(3, shinyjs::disabled(numericInput(inputId=ns("s_param2"), label="S",   value=NA_real_, step=0.1))),
          ),
          fluidRow(
            column(6, mod_reference_ui(id=ns("reference"))),
            column(6, numericInput(inputId=ns("l_param"), label="L", value=3, step=1)),
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p1"), label="p1", selected=1e-4, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("p2"), label="p2", selected=1e-4, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("p12"), label="p12", selected=1e-5, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("h4"), label="H4 decision", selected=0.5, choices=seq(0,1,by=0.1), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("r2thr"), label="R2 thresh", selected=0.01, choices=c(0.0001,0.001,0.01,seq(0.1,0.9,0.1),0.9999), grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("p_secondary"), label="p_secondary", selected=1e-6, choices=c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0), grid=TRUE))
          ),
          fluidRow(
            column(6, sliderTextInput(inputId=ns("maf_thresh"), label="EAF thresh", choices= c(0.0,0.01,0.05,0.1,0.2), selected=0.0, grid=TRUE)),
            column(6, sliderTextInput(inputId=ns("coverage"), label="Coverage %", selected=0.95, choices=c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0), grid=TRUE))
          )
        ) # end tagList

        reference_mod <<- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=TRUE)

      }

      return(controls)
    })


    #==========================================
    # function / method choice and input control activation
    #==========================================
    session$userData[[ns("source_type1")]] <- observeEvent(input$source_type1, {

      shinyjs::toggleState(id="s_param1", condition = !input$source_type1=="quant")
      shinyjs::toggleState(id="sd_y1",    condition =  input$source_type1=="quant")

    })
    session$userData[[ns("source_type2")]] <- observeEvent(input$source_type2, {

      shinyjs::toggleState(id="s_param2", condition = !input$source_type2=="quant")
      shinyjs::toggleState(id="sd_y2",    condition =  input$source_type2=="quant")

    })
    session$userData[[ns("sample_size1")]] <- observeEvent(list(source_1_module$data, input$n_param1), {

      req(source_1_module$data, !is.null(input$n_param1), !is.null(input$sd_y1))

      # get sample size
      if("N" %in% names(source_1_module$data) && is.na(input$n_param1)) {
        sample_size <- max(source_1_module$data$N, na.rm=TRUE)
        updateNumericInput(inputId="n_param1", value=sample_size)
      } else {
        sample_size <- input$n_param1
      }

      # calculate sdY if not given
      if(is.na(input$sd_y1) && all(c("SE","EAF") %in% names(source_1_module$data)) && !is.na(sample_size)) {
        if(input$source_type1 == "quant") {
          sdY_estimate <- sdY.est(source_1_module$data$SE ^ 2,
                                  source_1_module$data$EAF,
                                  sample_size)
          updateNumericInput(inputId="sd_y1", value=sdY_estimate)
        } else {
          updateNumericInput(inputId="sd_y1", value=NA_real_)
        }
      }

    })
    session$userData[[ns("sample_size2")]] <- observeEvent(list(source_2_module$data, input$n_param2), {

      req(source_2_module$data, !is.null(input$n_param2), !is.null(input$sd_y2))

      # get sample size
      if("N" %in% names(source_2_module$data) && is.na(input$n_param2)) {
        sample_size <- max(source_2_module$data$N, na.rm=TRUE)
        updateNumericInput(inputId="n_param2", value=sample_size)
      } else {
        sample_size <- input$n_param2
      }

      # calculate sdY if not given
      if(is.na(input$sd_y2) && all(c("SE","EAF") %in% names(source_2_module$data)) && !is.na(sample_size)) {
        if(input$source_type2 == "quant") {
          sdY_estimate <- sdY.est(source_2_module$data$SE ^ 2,
                                  source_2_module$data$EAF,
                                  sample_size)
          updateNumericInput(inputId="sd_y2", value=sdY_estimate)
        } else {
          updateNumericInput(inputId="sd_y2", value=NA_real_)
        }
      }

    })
    session$userData[[ns("method")]] <- observeEvent(list(input$method), {

      req(input$method)

      if(input$method == "mask") {
        updatePrettyRadioButtons(inputId="mode", selected="iterative")
      }
      shinyjs::toggleState(id="l_param",     condition = !input$method == "single")
      shinyjs::toggleState(id="r2thr",       condition =  input$method == "mask")
      shinyjs::toggleState(id="mode",        condition = !input$method == "mask")
      shinyjs::toggleState(id="p_secondary", condition = !input$method == "single")
      shinyjs::toggleState(id="p_secondary", condition = !input$method == "single")
      shinyjs::toggleState(id="reference-reference", condition= input$method=="cond")

    })



  })
}



############## Helper functions ################



#' @title Coloc package standard dataset
#' @description
#' Create the standard dataset for the [coloc] package
#' @param dat a data.table, with at least columns c("RSID","BP","BETA","SE") +/- c("EAF","N"), #"EA","OA",
#' @param type a string, either c("quant","cc")
#' @param sdY a numeric, std of the variable
#' @param N an integer, the sample size
#' @param S a numeric, proportion of cases
#' @param ld an LD matrix
#' @return a coloc dataset
#' @export
#'
make_coloc_dataset <- function(dat, type, sdY=NA_real_, N=NA_real_, S=NA_real_, ld=NULL) {

  # copy
  dat <- data.table::copy(dat)

  # checks
  type <- match.arg(type, choices = c("quant","cc"))
  stopifnot("Standard names missing from the input to make a coloc dataset" = all(c("RSID","BP","BETA","SE") %in% names(dat))) #"EA","OA",
  if(type=="quant") {
    if(is.na(sdY) && !(all(c("EAF","N") %in% names(dat)) | all(c(!is.na(N), "EAF" %in% names(dat))))) {
      showNotification("`sdY` or data columns `EAF` and `N` (or input value `N`) must be provided", type="error")
      return(NULL)
    }
  }

  # if using SuSiE LD based method, use alleles for RSID
  if(!is.null(ld)) {

    if(!"LD_RSID_allele" %in% names(dat)) {
      showNotification("If running SuSiE finemapping please provide `dat` and `ld` output from genepi.utils::ld_matrix()", type="error")
      return(NULL)
    }

    # make multiallelic variants unique ids
    dat[, RSID_local := LD_RSID_allele]
  } else {

    dat[, RSID_local := RSID]

  }

  # test for and deal with duplicate IDs
  if(sum(duplicated(dat$RSID_local), na.rm=TRUE) > 0) {
    dat <- dat[dat[, .I[which.min(P)], by=RSID_local]$V1]
    showNotification("Duplicate RSIDs found in source, taking the lowest P-value variant", type="warning")
  }

  # test for and deal with NAs
  if(sum(rowSums(is.na(dat[,list(RSID_local,BP,BETA,SE,EAF)]))) > 0) {

    na_rows <- rowSums(is.na(dat[,list(RSID_local,BP,BETA,SE,EAF)])) > 0
    dat <- dat[!na_rows, ]
    showNotification(paste0(sum(na_rows), " rows with NA values removed"), type="warning")

  }

  # create the base coloc object
  D1 <- list(
    snp      = dat$RSID_local,
    position = dat$BP,
    beta     = dat$BETA,
    varbeta  = dat$SE ^ 2,
    type     = type
  )

  # add additional parameters if provided
  if("EAF" %in% names(dat) && !is.na(N)) {

    D1$MAF <- dat$EAF
    D1$N   <- N

  } else if (all(c("EAF","N") %in% names(dat))) {

    D1$MAF <- dat$EAF
    D1$N   <- max(dat$N, na.rm=T)

  } else {

    showNotification("`EAF` must be provided; also, N` must be provided either as a parameter, or as a column in `dat`", type="error")
    return(NULL)

  }

  # proportion of sample that are cases
  if(!is.na(S)) {

    D1$s <- S

  }

  # add the st.dev of the Y variable
  if(!is.na(sdY)) {

    D1$sdY <- sdY

  }

  # add the LD matrix
  if(!is.null(ld)) {

    D1$LD <- ld[D1$snp, D1$snp]

  }

  # check the dataset
  tryCatch({
    if(!is.null(ld)) {
      coloc::check_dataset(D1, req="LD")
    } else {
      coloc::check_dataset(D1)
    }
  },
  error=function(e) {
    showNotification(paste0("Coloc dataset check failed - ", e), type="error", duration=10)
    return(NULL)
  })

  # return the object
  return(D1)
}


#' @title Estimate trait variance, internal function
#' @description
#' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
#'
#' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
#' var(X) = 2*maf*(1-maf)
#' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
#'
#' @param vbeta vector of variance of coefficients
#' @param maf vector of MAF (same length as vbeta)
#' @param n sample size
#' @return estimated standard deviation of Y
#' @references https://github.com/chr1swallace/coloc/blob/a44c817137daa0f4d0de56c9b77693f8228ab5fc/R/claudia.R#L135C1-L157C2
#' @author Chris Wallace
#' @noRd
sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}



#' @title get_harmonised_ld
#' @param dat .
#' @param ref_path .
#' @param chr .
#' @param from .
#' @param to .
#' @param method .
#' @return .
#' @export
#'
get_harmonised_ld <- function(dat, ref_path, chr, from, to, method, colmap=NULL) {

  # get the reference file
  if(grepl("1kGv3", ref_path)) {

    plink_ref <- make_ref_subset(ref_path=ref_path, chrom=chr, from=from, to=to)
    plink2    <- get_plink2_exe()
    ukbb_ref  <- NULL

  } else if(grepl("UKB_LD", ref_path)) {

    # get the UKBB LD file path (downloads if not in cache)
    ukbb_ref_dt <- genepi.utils::download_ukbb_ld(chr=chr, bp_start=from, bp_end=to, ukbb_ld_cache=file.path(ref_path, "cache"))
    ukbb_ref    <- ukbb_ref_dt$root_file
    plink_ref   <- NULL
    plink2      <- NULL

  }

  # get LD matrix and harmonise the data
  dat_ld_obj <- genepi.utils::ld_matrix(dat=dat, colmap=colmap, method=method, plink2=plink2, plink_ref=plink_ref, ukbb_ref=ukbb_ref)

  return(dat_ld_obj)

}



coloc_to_plotting <- function(result, h4=0.5, coverage=0.95) {

  # coloc.abf
  if(!is.data.frame(result$summary)) {

    summary <- data.table::data.table(RSID_1 = result$results$snp[which.max(result$results$SNP.PP.H4)],
                                      RSID_2 = result$results$snp[which.max(result$results$SNP.PP.H4)],
                                      PPH4   = result$summary[['PP.H4.abf']],
                                      step   = factor(1))

    # coloc.signals & coloc.susie
  } else {

    summary <- data.table::data.table(RSID_1 = result$summary$hit1,
                                      RSID_2 = result$summary$hit2,
                                      PPH4   = result$summary$PP.H4.abf,
                                      step   = factor(1:nrow(result$summary)))

  }

  # filter hits by H4 threshold
  summary <- summary[PPH4 >= h4, ]
  summary[, uid := apply(summary[, list(RSID_1, RSID_2)], 1, function(x) paste(sort(x), collapse = "_"))]
  summary <- unique(summary, by="uid")

  if(nrow(summary)==0) {
    summary[, group := factor(0)]
    summary[, index := logical(0)] # empty
  } else {
    summary[, group := factor(step, labels=1:length(unique(step)))]
    summary[, index := TRUE]
  }

  return(summary[, list(RSID_1,RSID_2,group,index)])

  # not sure I know what the PP is referring to vs. what the functions say are the top hits.
  # some conditioning that isn't represented in the PP values for each step??
  # Therefore, just plot what it says is the top hits, rather than try work out which points are in the
  # credible sets.
  #
  # # results
  # results <- data.table::as.data.table(result$results) |>
  #   data.table::melt(id.vars       = "snp",
  #                    measure.vars  = grep("SNP.PP", names(result$results), value=TRUE),
  #                    variable.name = "step",
  #                    value.name    = "PP")
  # data.table::setnames(results, "snp", "RSID")
  # data.table::setorder(results, step, -PP)
  # results[, step := factor(step, labels=1:length(unique(step)))]
  # results <- results[step %in% summary$step, ]
  # results[, CUMSUM_PP := cumsum(PP), by="step"]
  # results[ , CS := ifelse(seq_len(.N) <= which(CUMSUM_PP >= coverage)[1], TRUE, FALSE), by="step"]
  # results <- results[CS==TRUE, ]
  # results[results[, .I[which.min(CUMSUM_PP)], by=step]$V1, index := TRUE]
  # if(nrow(results)==0) {
  #   results[, group := factor(step)]
  # } else {
  #   results[, group := factor(step, labels=1:length(unique(step)))]
  # }
  #
  # # join
  # results[summary, c('RSID_1','RSID_2') := list(RSID_1,RSID_2), on="step"]
  #
  # return(results[, list(RSID,RSID_1,RSID_2,group,index)])

}
