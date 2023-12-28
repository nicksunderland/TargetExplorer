#' step UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets datamods
#'
mod_coloc_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_coloc_ui(ns=",ns("foo"),")"))
  # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
  div(
    id = id,
    sidebarLayout(position = "right",
                  sidebarPanel(width = 3,
                               fluidRow(column(10, p(strong(paste0("Controls [",id,"]")))),
                                        column(2,  mod_remove_ui(id=ns("remove")))),
                               hr(),
                               fluidRow(
                                 column(6,
                                        mod_source_select_ui(ns("source_1")),
                                        fluidRow(
                                          column(6,
                                                 prettyRadioButtons(inputId = ns("source_1_type"),
                                                                    label   = "Type",
                                                                    choices = c("quant","cc"),
                                                                    selected = "quant")),
                                          column(6,
                                                 numericInput(inputId = ns("sd_y1"),
                                                              label = "sdY",
                                                              value = NULL,
                                                              step  = 0.1))
                                        )),
                                 column(6,
                                        mod_source_select_ui(ns("source_2")),
                                        fluidRow(
                                          column(6,
                                                 prettyRadioButtons(inputId = ns("source_2_type"),
                                                                    label   = "Type",
                                                                    choices = c("quant","cc"),
                                                                    selected = "quant")),
                                          column(6,
                                                 numericInput(inputId = ns("sd_y2"),
                                                              label = "sdY",
                                                              value = NULL,
                                                              step  = 0.1))
                                        ))
                               ),
                               hr(),
                               fluidRow(
                                 column(6,
                                        sliderTextInput(inputId  = ns("coloc_p1"),
                                                        label    = "p1",
                                                        choices  = c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0),
                                                        selected = 1e-4,
                                                        grid     = TRUE),
                                        sliderTextInput(inputId  = ns("coloc_p12"),
                                                        label    = "p12",
                                                        choices  = c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0),
                                                        selected = 1e-5,
                                                        grid     = TRUE)),
                                 column(6,
                                        sliderTextInput(inputId  = ns("coloc_p2"),
                                                        label    = "p2",
                                                        choices  = c(5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1.0),
                                                        selected = 1e-4,
                                                        grid     = TRUE),
                                        sliderTextInput(inputId  = ns("coloc_h4"),
                                                        label    = "H4 decision",
                                                        choices  = seq(0,1,by=0.1),
                                                        selected = 0.5,
                                                        grid     = TRUE),
                                 )
                               ),
                               fluidRow(
                                 column(6,
                                        mod_reference_ui(id=ns("reference")),
                                        actionButton(inputId = ns("run"),
                                                     label   = "Run")),
                                 column(6,
                                        prettyRadioButtons(inputId  = ns("assumption"),
                                                           label    = "Assumption",
                                                           choices  = c("Single","SuSiE"),
                                                           selected = "Single",
                                                           inline   = TRUE),
                                        prettyRadioButtons(inputId  = ns("method"),
                                                           label    = "Method",
                                                           choices  = c("Finemap","Coloc"),
                                                           selected = "Finemap",
                                                           inline   = TRUE))
                               ),
                               hr(),
                               fluidRow(
                                 column(6,
                                        sliderTextInput(inputId  = ns("credible_set_p"),
                                                        label    = "Credible set %",
                                                        choices  = c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0),
                                                        selected = 0.95,
                                                        grid     = TRUE)),
                                 column(6,
                                        selectInput(inputId = ns("downstream_dataset"),
                                                    label   = "Downstream dataset",
                                                    choices = c("Dataset 1", "Dataset 2")))
                               )
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            p(strong("Colocalisation:")),
                            fluidRow(
                              column(6,
                                     plotOutput(outputId = ns("locus_plot"), height = "500px"),
                              ),
                              column(6,
                                     plotOutput(outputId = ns("prob_plot1"), height = "250px"),
                                     plotOutput(outputId = ns("prob_plot2"), height = "250px"),
                                     tableOutput(outputId = ns("coloc_plot_table"))
                              )
                            ),

                  ) # main panel end
    ), # sidebar layout end
    hr()
  ) # div end
}


#' coloc Server Functions
#' @import ggplot2 ggrepel
#' @importFrom coloc coloc.abf
#' @noRd
mod_coloc_server <- function(id, app){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_coloc_server(ns=",ns("foo"),")"))

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- nlog10P <- position <- tissue_label <- SNP.PP.H4 <- P <- NULL
    TISSUE <- index <- NULL

    #==========================================
    # Data, remove, and reference module servers for the Coloc module
    #==========================================
    data_mod      <- mod_data_server(id="data", gene_module=app$modules$gene)
    remove_mod    <- mod_remove_server(id="remove", app=app, parent_id=id, parent_inputs=input)
    reference_mod <- mod_reference_server(id="reference", label="LD reference", gene_module=app$modules$gene)
    source_1_mod  <- mod_source_select_server(id="source_1", app=app, source_type=c("GWAS","eQTL"), label="Dataset 1")
    source_2_mod  <- mod_source_select_server(id="source_2", app=app, source_type=c("GWAS","eQTL"), label="Dataset 2")


    #==========================================
    # Run finemap button
    #==========================================
    session$userData[[ns("run")]] <- observeEvent(input$run, {

      # catch coloc package warnings and print as notifications
      tryCatch({

        # DATASET 1 FINEMAP
        if(input$method=="Finemap" && input$assumption == "Single" && input$downstream_dataset == "Dataset 1" && !is.null(source_1_mod$data)) {

          # make and check the dataset
          D1 <- make_coloc_dataset(dat  = source_1_mod$data,
                                   type = input$source_1_type,
                                   sdY  = input$sd_y1,
                                   ld   = NULL)

          # run the coloc finemap function
          data_mod$data2$results <- coloc::finemap.abf(D1, p1=input$coloc_p1)
          data_mod$source <- "Dataset 1"
          data_mod$data  <- calc_credible_set(rsids          = data_mod$data2$results$snp,
                                              snp.pp         = data_mod$data2$results$SNP.PP,
                                              dat_join_to    = data.table::copy(source_1_mod$data),
                                              credible_set_p = input$credible_set_p)

          # DATASET 1 SuSiE - finemap
        } else if(input$method=="Finemap" && input$assumption == "SuSiE" && input$downstream_dataset == "Dataset 1" && !is.null(source_1_mod$data)) {

          # get the reference file
          plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                       chrom    = app$modules$gene$chr,
                                       from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
                                       to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)

          # create the LD matrix for the region variants
          dat_ld_obj <- genepi.utils::ld_matrix(variants     = source_1_mod$data,
                                                with_alleles = TRUE,
                                                plink2       = get_plink2_exe(),
                                                plink_ref    = plink_ref)

          # make and check the dataset 1
          D1 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
                                   type = input$source_1_type,
                                   sdY  = input$sd_y1,
                                   ld   = dat_ld_obj[["ld_mat"]])

          # run SuSie
          data_mod$data2  <- coloc::runsusie(D1, coverage=input$credible_set_p)
          data_mod$data2$kriging_rss <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), dat_ld_obj[["ld_mat"]], n=D1$N)$conditional_dist
          data_mod$data2$kriging_rss$RSID <- D1$snp
          data_mod$data2$kriging_rss$outlier <- ifelse(data_mod$data2$kriging_rss$logLR  > 2 &
                                                         abs(data_mod$data2$kriging_rss$z) > 2, TRUE, FALSE)
          data_mod$source <- "Dataset 1"
          data_mod$data   <- calc_credible_set(rsids          = names(data_mod$data2$pip),
                                               snp.pp         = data_mod$data2$pip,
                                               dat_join_to    = data.table::copy(source_1_mod$data),
                                               credible_sets  = data_mod$data2$sets$cs,
                                               credible_set_p = NULL)

          # DATASET 2 FINEMAP
        } else if(input$method=="Finemap" && input$assumption == "Single" && input$downstream_dataset == "Dataset 2" && !is.null(source_2_mod$data)) {

          # make and check the dataset
          D2 <- make_coloc_dataset(dat  = source_2_mod$data,
                                   type = input$source_2_type,
                                   sdY  = input$sd_y2,
                                   ld   = NULL)

          # run the coloc finemap function
          data_mod$data2$results <- coloc::finemap.abf(D2, p1=input$coloc_p1)
          data_mod$source <- "Dataset 2"
          data_mod$data  <- calc_credible_set(rsids          = data_mod$data2$results$snp,
                                              snp.pp         = data_mod$data2$results$SNP.PP,
                                              dat_join_to    = data.table::copy(source_1_mod$data),
                                              credible_set_p = input$credible_set_p)

          # DATASET 2 SuSiE - finemap
        } else if(input$method=="Finemap" && input$assumption == "SuSiE" && input$downstream_dataset == "Dataset 2" && !is.null(source_2_mod$data)) {

          # get the reference file
          plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                       chrom    = app$modules$gene$chr,
                                       from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
                                       to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)

          # create the LD matrix for the region variants
          dat_ld_obj <- genepi.utils::ld_matrix(variants     = source_2_mod$data,
                                                with_alleles = TRUE,
                                                plink2       = get_plink2_exe(),
                                                plink_ref    = plink_ref)

          # make and check the dataset 1
          D2 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
                                   type = input$source_2_type,
                                   sdY  = input$sd_y2,
                                   ld   = dat_ld_obj[["ld_mat"]])

          # run SuSie
          data_mod$data2  <- coloc::runsusie(D2, coverage=input$credible_set_p)
          data_mod$data2$kriging_rss <- susieR::kriging_rss(D2$beta/(D2$varbeta ^ 0.5), dat_ld_obj[["ld_mat"]], n=D2$N)$conditional_dist
          data_mod$data2$kriging_rss$RSID <- D2$snp
          data_mod$data2$kriging_rss$outlier <- ifelse(data_mod$data2$kriging_rss$logLR  > 2 &
                                                         abs(data_mod$data2$kriging_rss$z) > 2, TRUE, FALSE)
          data_mod$source <- "Dataset 2"
          data_mod$data   <- calc_credible_set(rsids          = names(data_mod$data2$pip),
                                               snp.pp         = data_mod$data2$pip,
                                               dat_join_to    = data.table::copy(source_2_mod$data),
                                               credible_sets  = data_mod$data2$sets$cs,
                                               credible_set_p = NULL)


          # DATASET 1&2 COLOC
        } else if(input$method=="Coloc" && input$assumption == "Single" && !is.null(source_1_mod$data) && !is.null(source_2_mod$data)) {

          # harmonise the datasets
          harm <- genepi.utils::harmonise(gwas1 = source_1_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
                                          gwas2 = source_2_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
                                          gwas1_trait = "d1",
                                          gwas2_trait = "d2",
                                          merge = c("SNP"="SNP"))

          # make and check the dataset 1
          D1 <- make_coloc_dataset(dat  = harm[, list(RSID=RSID_d1,BP=BP_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
                                   type = input$source_1_type,
                                   sdY  = input$sd_y1,
                                   ld   = NULL)

          # make and check the dataset 2
          D2 <- make_coloc_dataset(dat  = harm[, list(RSID=RSID_d2,BP=BP_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
                                   type = input$source_2_type,
                                   sdY  = input$sd_y2,
                                   ld   = NULL)

          # run the coloc abf function
          res <- coloc::coloc.abf(D1, D2,
                                  p1  = input$coloc_p1,
                                  p2  = input$coloc_p2,
                                  p12 = input$coloc_p12)

          # assign the results
          data_mod$data2 <- res
          data_mod$source <- "Coloc"

          # populate data_mod$data with whichever dataset we plan on feeding downstream
          if(input$downstream_dataset == "Dataset 1") {

            dat_join_to = source_1_mod$data

          } else if(input$downstream_dataset == "Dataset 2") {

            dat_join_to = source_2_mod$data

          }
          data_mod$data <- calc_credible_set(rsids          = data_mod$data2$results$snp,
                                             snp.pp         = data_mod$data2$results$SNP.PP.H4,
                                             dat_join_to    = dat_join_to,
                                             credible_set_p = input$credible_set_p)

          # DATASET 1&2 Coloc SuSiE
        } else if(input$method=="Coloc" && input$assumption == "SuSiE" && !is.null(source_1_mod$data) && !is.null(source_2_mod$data)) {

          # harmonise the datasets
          harm <- genepi.utils::harmonise(gwas1 = source_1_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
                                          gwas2 = source_2_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
                                          gwas1_trait = "d1",
                                          gwas2_trait = "d2",
                                          merge = c("SNP"="SNP"))

          # get the reference file
          plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                       chrom    = app$modules$gene$chr,
                                       from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
                                       to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)

          # create the LD matrix for the region variants
          dat_ld_obj1 <- genepi.utils::ld_matrix(variants     = harm[, list(RSID=RSID_d1,BP=BP_d1,EA=EA_d1,OA=OA_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
                                                 with_alleles = TRUE,
                                                 plink2       = get_plink2_exe(),
                                                 plink_ref    = plink_ref)
          dat_ld_obj2 <- genepi.utils::ld_matrix(variants     = harm[, list(RSID=RSID_d2,BP=BP_d2,EA=EA_d2,OA=OA_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
                                                 with_alleles = TRUE,
                                                 plink2       = get_plink2_exe(),
                                                 plink_ref    = plink_ref)

          # make and check the dataset 1
          D1 <- make_coloc_dataset(dat  = dat_ld_obj1[["dat"]],
                                   type = input$source_1_type,
                                   sdY  = input$sd_y1,
                                   ld   = dat_ld_obj1[["ld_mat"]])

          # make and check the dataset 2
          D2 <- make_coloc_dataset(dat  = dat_ld_obj2[["dat"]],
                                   type = input$source_2_type,
                                   sdY  = input$sd_y2,
                                   ld   = dat_ld_obj2[["ld_mat"]])


          # run SuSiE
          D1_susie <- coloc::runsusie(D1, coverage=input$credible_set_p)
          D2_susie <- coloc::runsusie(D2, coverage=input$credible_set_p)

          # run the coloc abf function
          res <- coloc::coloc.susie(D1_susie, D2_susie)


          # Browse[1]>  print(res$summary)
          # nsnps            hit1            hit2     PP.H0.abf     PP.H1.abf     PP.H2.abf    PP.H3.abf     PP.H4.abf  idx1  idx2
          # <int>          <char>          <char>         <num>         <num>         <num>        <num>         <num> <int> <int>
          #   1:   390  rs10305435_C_T  rs10305435_C_T 2.190624e-175  2.093143e-89  2.093143e-89 0.000000e+00  1.000000e+00     1     1
          # 2:   390  rs10305445_C_T  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  1.628576e-88     2     1
          # 3:   390   rs2300613_G_A  rs10305435_C_T 5.057475e-161  1.046572e-86  4.832421e-75 1.000000e+00  9.458692e-46     3     1
          # 4:   390   rs2300614_A_G  rs10305435_C_T 6.092174e-263  1.046572e-86 5.821077e-177 1.000000e+00  8.832468e-88     4     1
          # 5:   390    rs877446_A_G  rs10305435_C_T 3.379615e-298  1.046572e-86 3.229225e-212 1.000000e+00  1.935952e-85     6     1
          # 6:   390  rs10305442_A_G  rs10305435_C_T 2.369280e-277  1.046572e-86 2.263849e-191 1.000000e+00  1.214910e-88     7     1
          # 7:   390 rs115493670_C_T  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  4.568936e-89     8     1
          # 8:   390   rs9283907_A_G  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  1.434933e-88     9     1
          # 9:   390   rs9283906_T_C  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  4.752047e-89    10     1
          # 10:   390   rs2235868_A_C  rs10305435_C_T 3.821291e-149  1.046572e-86  3.651246e-63 1.000000e+00  3.825684e-25     5     1
          # 11:   390  rs10305435_C_T  rs10305445_C_T  0.000000e+00  0.000000e+00  1.046572e-86 1.000000e+00  1.628576e-88     1     2

          #
          #
          # # assign the results
          # data_mod$data2 <- res
          # data_mod$source <- "Coloc"
          #
          # # populate data_mod$data with whichever dataset we plan on feeding downstream
          # if(input$downstream_dataset == "Dataset 1") {
          #
          #   dat_join_to = source_1_mod$data
          #
          # } else if(input$downstream_dataset == "Dataset 2") {
          #
          #   dat_join_to = source_2_mod$data
          #
          # }
          # data_mod$data <- calc_credible_set(rsids          = data_mod$data2$results$snp,
          #                                    snp.pp         = data_mod$data2$results$SNP.PP.H4,
          #                                    dat_join_to    = dat_join_to,
          #                                    credible_set_p = input$credible_set_p)
          #

        }

        # button text back to black
        shinyjs::runjs(paste0('document.getElementById("', ns("run"), '").style.color = "black";'))

      },
      error=function(e) {
        showNotification(paste0("coloc::finemap.abf() failed - ", e), type="error", duration=10)
        return(NULL)
      })


    })


    #==========================================
    # Observe controls - make run button red if changed
    #==========================================
    observeEvent(list(input$downstream_dataset,
                      input$source_1,
                      input$source_2,
                      input$source_1_type,
                      input$source_2_type,
                      input$sd_y1,
                      input$sd_y2,
                      input$coloc_p1,
                      input$coloc_p12,
                      input$coloc_p2,
                      input$coloc_h4,
                      input$ld_reference,
                      input$assumption,
                      input$method,
                      input$credible_set_p
    ), {

      shinyjs::runjs(paste0('document.getElementById("', ns("run"), '").style.color = "red";'))

    }, ignoreInit=TRUE)


    #==========================================
    # Observe credible set P-value threshold
    #==========================================
    session$userData[[ns("credible_set_p")]] <- observeEvent(input$credible_set_p, {

      if(is.null(data_mod$data)) return(NULL)
      req("cumsum_pp_h4" %in% names(data_mod$data))

      # determine credible set - flag as coloc
      data_mod$data[, coloc := ifelse(cumsum_pp_h4 > input$credible_set_p, TRUE, FALSE)]

    })


    #==========================================
    # Observe additions / deletions of modules
    #==========================================
    session$userData[[ns("app-modules")]] <- observeEvent(app$modules, {
      updateSelectInput(session, inputId="source_1", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
      updateSelectInput(session, inputId="source_2", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
    })


    #==========================================
    # Observe input selection for choices init as observing just app$modules doesn't work at module init...
    #==========================================
    session$userData[[ns("sources")]] <- observeEvent(list(input$source_1,input$source_2), {
      if(is.null(input$source_1) || input$source_1=="") {
        updateSelectInput(session, inputId="source_1", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
      }
      if(is.null(input$source_2) || input$source_2=="") {
        updateSelectInput(session, inputId="source_2", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
      }
    })


    #==========================================
    # Observe method input (SuSiE / single causal variant assumption)
    #==========================================
    session$userData[[ns("assumption")]] <- observeEvent(input$assumption, {
      if(input$assumption == "SuSiE") {
        shinyjs::enable("ld_reference")
      } else if(input$assumption == "Single") {
        shinyjs::disable("ld_reference")
      }
    })


    #==========================================
    # Results table
    #==========================================
    output$coloc_plot_table <- renderTable({
      req(data_mod$data2$summary)
      return(t(as.data.frame(data_mod$data2$summary)))
    })


    #==========================================
    # Locus plot
    #==========================================
    output$locus_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(source_1_mod$data) | !is.null(source_2_mod$data), 'No data imported for either dataset')
      )

      # colours
      color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

      # base plot
      p <- ggplot() +
        theme_classic() +
        lims(x = c(app$modules$gene$start - app$modules$gene$flanks_kb*1000, app$modules$gene$end + app$modules$gene$flanks_kb*1000)) +
        labs(x = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y = expression(paste("-log"[10], plain(P)))) +
        theme(legend.position="top")

      # y axis max and min
      y_max <- c(NA_real_, NA_real_)

      # plot source 1 data positive
      if(!is.null(source_1_mod$data)) {

        y_max[[1]] <- max(source_1_mod$data$nlog10P, na.rm=TRUE)

        p <- p +
          annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05)

        # add downstream label if UI indicates
        if(input$downstream_dataset=="Dataset 1") {
          p <- p +
            annotate(geom = "text", x=app$modules$gene$end+(app$modules$gene$flanks_kb*1000)-1e4, y=ceiling(y_max[[1]]), label="downstream", color="darkred", alpha=0.5)
        }

        # FINEMAPPING SOURCE 1 - plot the points
        if(!is.null(data_mod$data) && data_mod$source == "Dataset 1" && input$downstream_dataset == "Dataset 1") {

          p <- p +
            geom_point(data    = data_mod$data,
                       mapping = aes(x=BP, y=nlog10P, color=PP), shape=19) +
            geom_point(data    = data_mod$data[which(coloc==TRUE), ],
                       mapping = aes(x=BP, y=nlog10P, fill=credible_set), stroke=NA, shape=24, size=3) +
            viridis::scale_color_viridis(option="magma") +
            geom_label_repel(data    = data_mod$data[which(coloc==TRUE), ],
                             mapping = aes(label=RSID, x=BP, y=nlog10P),
                             max.overlaps = Inf) +
            labs(color="PP", fill=paste0(input$credible_set_p*100,"% credible sets"))

          # COLOC SOURCE 1 - plot the points
        } else if(!is.null(data_mod$data) && data_mod$source == "Coloc") {

          # if data_mod$source == "Coloc"
          # data_mod$data is either the upstream source 1 or source 2, with the joined coloc data - extra cols:
          # PP<num>, coloc<lgl>, credible_set<fct>

          # join to data source
          plot_data <- source_1_mod$data[data_mod$data[, list(RSID,PP,coloc,credible_set)], on="RSID"]
          plot_data[source_2_mod$data, nlog10P_lower := i.nlog10P, on="RSID"]

          p <- p +
            geom_point(data    = plot_data,
                       mapping = aes(x=BP, y=nlog10P, color=PP), shape=19) +
            geom_segment(data    = plot_data[which(coloc==TRUE), ],
                         mapping = aes(x=BP,xend=BP,y=nlog10P,yend=-nlog10P_lower), linetype="dotted", color="darkred") +
            geom_point(data    = plot_data[which(coloc==TRUE), ],
                       mapping = aes(x=BP, y=nlog10P, fill=credible_set), stroke=NA, shape=24, size=3) +
            viridis::scale_color_viridis(option="magma") +
            geom_label_repel(data    = plot_data[which(coloc==TRUE), ],
                             mapping = aes(label=RSID, x=BP, y=nlog10P),
                             max.overlaps = Inf) +
            labs(color="PP", fill=paste0(input$credible_set_p*100,"% credible sets"))


          # RAW SOURCE 1 - plot the points
        } else {

          p <- p +
            geom_point(data    = source_1_mod$data,
                       mapping = aes(x=BP, y=nlog10P),
                       color   = "lightgray")

        }

      }

      # plot source 2 data negative
      if(!is.null(source_2_mod$data)) {

        y_max[[2]] <- max(source_2_mod$data$nlog10P, na.rm=TRUE)

        p <- p +
          annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=-Inf, ymax=0, fill="blue", alpha = 0.05)

        # add downstream label if UI indicates
        if(input$downstream_dataset=="Dataset 2") {
          p <- p +
            annotate(geom = "text", x=app$modules$gene$end+(app$modules$gene$flanks_kb*1000)-1e4, y=-ceiling(y_max[[2]]), label="downstream", color="darkred", alpha=0.5)
        }

        # FINEMAPPING SOURCE 2 - plot the points
        if(!is.null(data_mod$data) && data_mod$source == "Dataset 2" &&  input$downstream_dataset == "Dataset 2") {

          p <- p +
            geom_point(data    = data_mod$data,
                       mapping = aes(x=BP, y=-nlog10P, color=PP), shape=19) +
            geom_point(data    = data_mod$data[which(coloc==TRUE), ],
                       mapping = aes(x=BP, y=-nlog10P, fill=credible_set), stroke=NA, shape=24, size=3) +
            viridis::scale_color_viridis(option="magma") +
            geom_label_repel(data    = data_mod$data[which(coloc==TRUE), ],
                             mapping = aes(label=RSID, x=BP, y=-nlog10P),
                             max.overlaps = Inf) +
            labs(color="PP", fill=paste0(input$credible_set_p*100,"% credible sets"))

          # COLOC SOURCE 2 - plot the points
        } else if(!is.null(data_mod$data) && data_mod$source == "Coloc") {

          # if data_mod$source == "Coloc"
          # data_mod$data is either the upstream source 1 or source 2, with the joined coloc data - extra cols:
          # PP<num>, coloc<lgl>, credible_set<fct>

          # join to data source
          plot_data <- source_2_mod$data[data_mod$data[, list(RSID,PP,coloc,credible_set)], on="RSID"]

          p <- p +
            geom_point(data    = plot_data,
                       mapping = aes(x=BP, y=-nlog10P, color=PP), shape=19) +
            geom_point(data    = plot_data[which(coloc==TRUE), ],
                       mapping = aes(x=BP, y=-nlog10P, fill=credible_set), stroke=NA, shape=24, size=3) +
            viridis::scale_color_viridis(option="magma") +
            labs(color="PP", fill=paste0(input$credible_set_p*100,"% credible sets"))


          # RAW SOURCE 2 - plot the points
        } else {

          p <- p +
            geom_point(data    = source_2_mod$data,
                       mapping = aes(x=BP, y=-nlog10P),
                       color   = "darkgray")

        }
      }

      # if two sets of data
      if(!is.null(source_1_mod$data) && !is.null(source_2_mod$data)) {

        # draw the x-axis at y=0
        p <- p + geom_hline(yintercept = 0)

      }

      # correct P values on the Y axis
      p <- p +
        scale_y_continuous(breaks =     seq(ifelse(is.na(y_max[[2]]), 0, -ceiling(y_max[[2]])),
                                            ifelse(is.na(y_max[[1]]), 0,  ceiling(y_max[[1]])), 1),
                           labels = abs(seq(ifelse(is.na(y_max[[2]]), 0, -ceiling(y_max[[2]])),
                                            ifelse(is.na(y_max[[1]]), 0,  ceiling(y_max[[1]])), 1)))


      return(p)
    })


    #==========================================
    # Prior probabilities plot
    #==========================================
    output$prob_plot1 <- renderPlot({

      # check data
      if(is.null(data_mod$data2$summary) & is.null(data_mod$data2$kriging_rss)) {

        p <- NULL

      } else if(!is.null(data_mod$data2$summary)) {

        p <- genepi.utils::plot_coloc_probabilities(data_mod$data2, rule=paste0("H4 > ", input$coloc_h4), type="prior")

      } else if(!is.null(data_mod$data2$kriging_rss)) {

        p <- ggplot(data    = data_mod$data2$kriging_rss,
                    mapping = aes(x=condmean, y=z)) +
          geom_abline(intercept = 0, slope = 1) +
          geom_point(mapping = aes(color=outlier)) +
          geom_label_repel(data = data_mod$data2$kriging_rss[data_mod$data2$kriging_rss$outlier==TRUE, ],
                           mapping = aes(label=RSID), max.overlaps = Inf) +
          scale_color_manual(values = c("TRUE"="red", "FALSE"="darkgrey")) +
          theme_classic() +
          labs(y = "Observed z scores", x = "Expected value") +
          theme(legend.position = "none")

      }

      return(p)
    })


    #==========================================
    # Posterior probabilities plot
    #==========================================
    output$prob_plot2 <- renderPlot({

      # check data
      validate(
        need(!is.null(data_mod$data2$summary), 'No colocation data found, run Coloc')
      )

      p <- genepi.utils::plot_coloc_probabilities(data_mod$data2, rule=paste0("H4 > ", input$coloc_h4), type="posterior")

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




make_coloc_dataset <- function(dat, type, sdY=NA_real_, ld=NULL) {

  # checks
  type <- match.arg(type, choices = c("quant","cc"))
  if(type=="quant") {
    stopifnot("Standard names missing from the input to make a coloc dataset" = all(c("RSID","BP","BETA","SE") %in% names(dat)))
    stopifnot("Standard data missing from the input to make a coloc dataset" = all(!is.null(sdY) | c("EAF","N") %in% names(dat)))
  } else {
    stopifnot("Standard names missing from the input to make a coloc dataset" = all(c("RSID","BP","BETA","SE") %in% names(dat)))
  }

  # make multiallelic variants unique ids
  if("RSID_allele" %in% names(dat)) {
    dat[, RSID := RSID_allele]
  }

  # test for and deal with duplicate IDs
  if(sum(duplicated(dat$RSID), na.rm=TRUE) > 0) {
    dat <- dat[dat[, .I[which.min(P)], by=RSID]$V1]
    showNotification("Duplicate RSIDs found in source 1, taking the lowest P-value variant", type="warning")
  }

  # create the base coloc object
  D1 <- list(
    snp      = dat$RSID,
    position = dat$BP,
    beta     = dat$BETA,
    varbeta  = dat$SE ^ 2,
    type     = type
  )

  # add additional parameters if provided
  if(all(c("EAF","N") %in% names(dat))) {
    D1$MAF <- dat$EAF
    D1$N   <- max(dat$N, na.rm=T)
  }
  if(!is.na(sdY)) {
    D1$sdY <- sdY
  }
  if(!is.null(ld)) {
    D1$LD <- ld
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




calc_credible_set <- function(rsids, snp.pp, dat_join_to, credible_sets=NULL, credible_set_p=NULL) {

  # checks
  stopifnot("Either pass credible_sets (coloc::runsusie [result]$sets$cs output) or credible_set_p, not both" = sum(sapply(list(credible_sets, credible_set_p), is.null))==1)

  # table the rsids and posterior probability
  # table the vectors
  dat_pp <- data.table::data.table(RSID = sub("^(rs[0-9]+).*","\\1", rsids),
                                   PP = snp.pp)

  # susie output
  if(!is.null(credible_sets)) {

    # extract the credible sets to dt
    cs <- data.table::data.table(RSID          = sub("L[0-9]+\\.(rs[0-9]+).*","\\1",names(unlist(credible_sets))),
                                 credible_set  = sub("L([0-9]+)\\..*","\\1",names(unlist(credible_sets))))

    # credible set factor
    cs[, credible_set := factor(credible_set, levels=sort(unique(as.integer(credible_set))))]

    # join to PP
    dat_pp[cs, credible_set := i.credible_set, on="RSID"]

    # add coloc flag
    dat_pp[, coloc := ifelse(!is.na(credible_set), TRUE, FALSE)]

    # standard signle variant assumption
  } else if(!is.null(credible_set_p)) {

    # order and get cumulative posterior probability
    dat_pp[order(PP, decreasing=TRUE), CUMSUM_PP := cumsum(PP)]

    # determine credible set - flag as coloc (i.e. when is the cumulative probability > 95% / whatever indicated)
    dat_pp[order(CUMSUM_PP), coloc := ifelse(.I <= which(CUMSUM_PP >= credible_set_p)[1], TRUE, FALSE)]

    # credible set factor
    dat_pp[, credible_set := ifelse(coloc, 1, NA)]
    dat_pp[, credible_set := factor(credible_set, levels=c(1))]

  }

  # join to the provided data
  dat_join_to[dat_pp, c("PP","coloc","credible_set") := list(i.PP, i.coloc, i.credible_set), on="RSID"]

  # return the data
  return(dat_join_to)
}




