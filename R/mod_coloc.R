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
  # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
  div(
    id = id,
    sidebarLayout(position = "right",
                  sidebarPanel(width = 3,
                               fluidRow(
                                 column(10, p(strong(paste0("Controls [",id,"]")))),
                                 column(2,  mod_remove_ui(id=ns("remove")))),
                               hr(),
                               fluidRow(
                                 column(6, mod_source_select_ui(ns("source_1"))),
                                 column(6, mod_source_select_ui(ns("source_2")))
                               ),
                               hr(),
                               fluidRow(
                                 column(4, p(strong("Colocalisation:"), style="text-align:right; margin-top: 6px;")),
                                 column(8,
                                        selectInput(inputId  = ns("package_select"),
                                                    label    = NULL,
                                                    choices  = c("Off","Common variants","r-coloc","r-susieR"),
                                                    selected = "Off"))
                               ),
                               uiOutput(ns("package_ui")),
                               fluidRow(
                                 column(4, p(strong("Output data:"), style="text-align:left; margin-top: 6px;")),
                                 column(8,
                                        selectInput(inputId = ns("downstream_dataset"),
                                                    label   = NULL,
                                                    choices = c("Dataset 1", "Dataset 2")))
                               )
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            fluidRow(
                              column(9, p(strong("Colocalisation:"))),
                              column(1, p(strong("Sensitivity:"), style="text-align:left; margin-top: 6px;")),
                              column(2, selectInput(inputId = ns("sensitivity_plot"), label=NULL, choices = c("Off","test"), selected = "Off"))
                            ),
                            # Locus zoom plot
                            uiOutput(outputId = ns("plot_area")),
                            tableOutput(outputId = ns("locus_plot_table")),
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


    #==========================================
    # Data, remove, and reference module servers for the Coloc module
    #==========================================
    data_mod      <- mod_data_server(id="data", gene_module=app$modules$gene)
    remove_mod    <- mod_remove_server(id="remove", app=app, parent_id=id, parent_inputs=input)
    source_1_mod  <- mod_source_select_server(id="source_1", app=app, source_type=c("GWAS"), label="Dataset 1")
    source_2_mod  <- mod_source_select_server(id="source_2", app=app, source_type=c("GWAS"), label="Dataset 2")


    #==========================================
    # Package select input - render package UI
    #==========================================
    output$package_ui <- renderUI({

      if(input$package_select == "Off") {

        controls <- fluidRow()
        data_mod$data  <- NULL
        data_mod$data2 <- NULL

      } else if(input$package_select == "Common variants") {

        controls <- fluidRow(p("TODO"))

      } else if(input$package_select == "r-coloc") {

        controls <- mod_r_coloc_ui(id=ns("r_coloc"))
        mod_r_coloc_server(id="r_coloc", gene_module=app$modules$gene, data_module=data_mod, functions=c("coloc.abf","coloc.signals","coloc.susie"))

      } else if(input$package_select == "r-susieR") {

        controls <- fluidRow(p("TODO"))
        # controls <- mod_r_susier_ui(id=ns("r_susier"))
        # mod_r_susier_server(id="r_susier", gene_module=app$modules$gene, data_module=data_mod)

      }

      # return the controls UI elements
      return(controls)
    })


    #==========================================
    # Plot area rendering - locus +/- sensitivity plots
    #==========================================
    output$plot_area <- renderUI({

      if(is.null(input$sensitivity_plot) || input$sensitivity_plot == "Off") {

        plot_area <- column(12, plotOutput(outputId=ns("locus_plot"), height="550px", brush=ns("locus_plot_brush")))

      } else {

        plot_area <- fluidRow(column(6, plotOutput(outputId=ns("locus_plot"), height="550px", brush=ns("locus_plot_brush"))),
                              column(6, plotOutput(outputId=ns("sensitivity_plot"), height="550px")))
      }
      return(plot_area)
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
    # Results table
    #==========================================
    output$locus_plot_table <- renderTable({

      req(data_mod$data2$summary)

      return(t(as.data.frame(data_mod$data2$summary)))

    })


    #==========================================
    # Sensitivity plots
    #==========================================
    output$sensitivity_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(input$sensitivity_plot) && input$sensitivity_plot != "Off", 'No sensitivity plots selected'),
        need(!is.null(data_mod$data), 'No colocalisation data found')
      )


      #--------------------------------------------
      # XXXXXXXXXX
      if(input$sensitivity_plot == "testing") {

        p <- NULL

      } else if(input$sensitivity_plot == "foo") {

        p <- genepi.utils::plot_coloc_probabilities(data_mod$data2, rule=paste0("H4 > ", input$coloc_h4), type="prior")

      } else if(input$sensitivity_plot == "bar") {

        p <- genepi.utils::plot_coloc_probabilities(data_mod$data2, rule=paste0("H4 > ", input$coloc_h4), type="posterior")

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














#
# #==========================================
# # Run finemap button
# #==========================================
# session$userData[[ns("run")]] <- observeEvent(input$run, {
#
#   # catch coloc package warnings and print as notifications
#   tryCatch({
#
#     # DATASET 1 FINEMAP
#     if(input$method=="Finemap" && input$assumption == "Single" && input$downstream_dataset == "Dataset 1" && !is.null(source_1_mod$data)) {
#
#       # make and check the dataset
#       D1 <- make_coloc_dataset(dat  = source_1_mod$data,
#                                type = input$source_1_type,
#                                sdY  = input$sd_y1,
#                                ld   = NULL)
#
#       # run the coloc finemap function
#       data_mod$data2$results <- coloc::finemap.abf(D1, p1=input$coloc_p1)
#       data_mod$source <- "Dataset 1"
#       data_mod$data  <- calc_credible_set(rsids          = data_mod$data2$results$snp,
#                                           snp.pp         = data_mod$data2$results$SNP.PP,
#                                           dat_join_to    = data.table::copy(source_1_mod$data),
#                                           credible_set_p = input$credible_set_p)
#
#       # DATASET 1 SuSiE - finemap
#     } else if(input$method=="Finemap" && input$assumption == "SuSiE" && input$downstream_dataset == "Dataset 1" && !is.null(source_1_mod$data)) {
#
#       # get the reference file
#       plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
#                                    chrom    = app$modules$gene$chr,
#                                    from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
#                                    to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)
#
#       # create the LD matrix for the region variants
#       dat_ld_obj <- genepi.utils::ld_matrix(variants     = source_1_mod$data,
#                                             with_alleles = TRUE,
#                                             plink2       = get_plink2_exe(),
#                                             plink_ref    = plink_ref)
#
#       # make and check the dataset 1
#       D1 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
#                                type = input$source_1_type,
#                                sdY  = input$sd_y1,
#                                ld   = dat_ld_obj[["ld_mat"]])
#
#       # run SuSie
#       data_mod$data2  <- coloc::runsusie(D1, coverage=input$credible_set_p)
#       data_mod$data2$kriging_rss <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), dat_ld_obj[["ld_mat"]], n=D1$N)$conditional_dist
#       data_mod$data2$kriging_rss$RSID <- D1$snp
#       data_mod$data2$kriging_rss$outlier <- ifelse(data_mod$data2$kriging_rss$logLR  > 2 &
#                                                      abs(data_mod$data2$kriging_rss$z) > 2, TRUE, FALSE)
#       data_mod$source <- "Dataset 1"
#       data_mod$data   <- calc_credible_set(rsids          = names(data_mod$data2$pip),
#                                            snp.pp         = data_mod$data2$pip,
#                                            dat_join_to    = data.table::copy(source_1_mod$data),
#                                            credible_sets  = data_mod$data2$sets$cs,
#                                            credible_set_p = NULL)
#
#       # DATASET 2 FINEMAP
#     } else if(input$method=="Finemap" && input$assumption == "Single" && input$downstream_dataset == "Dataset 2" && !is.null(source_2_mod$data)) {
#
#       # make and check the dataset
#       D2 <- make_coloc_dataset(dat  = source_2_mod$data,
#                                type = input$source_2_type,
#                                sdY  = input$sd_y2,
#                                ld   = NULL)
#
#       # run the coloc finemap function
#       data_mod$data2$results <- coloc::finemap.abf(D2, p1=input$coloc_p1)
#       data_mod$source <- "Dataset 2"
#       data_mod$data  <- calc_credible_set(rsids          = data_mod$data2$results$snp,
#                                           snp.pp         = data_mod$data2$results$SNP.PP,
#                                           dat_join_to    = data.table::copy(source_1_mod$data),
#                                           credible_set_p = input$credible_set_p)
#
#       # DATASET 2 SuSiE - finemap
#     } else if(input$method=="Finemap" && input$assumption == "SuSiE" && input$downstream_dataset == "Dataset 2" && !is.null(source_2_mod$data)) {
#
#       # get the reference file
#       plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
#                                    chrom    = app$modules$gene$chr,
#                                    from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
#                                    to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)
#
#       # create the LD matrix for the region variants
#       dat_ld_obj <- genepi.utils::ld_matrix(variants     = source_2_mod$data,
#                                             with_alleles = TRUE,
#                                             plink2       = get_plink2_exe(),
#                                             plink_ref    = plink_ref)
#
#       # make and check the dataset 1
#       D2 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
#                                type = input$source_2_type,
#                                sdY  = input$sd_y2,
#                                ld   = dat_ld_obj[["ld_mat"]])
#
#       # run SuSie
#       data_mod$data2  <- coloc::runsusie(D2, coverage=input$credible_set_p)
#       data_mod$data2$kriging_rss <- susieR::kriging_rss(D2$beta/(D2$varbeta ^ 0.5), dat_ld_obj[["ld_mat"]], n=D2$N)$conditional_dist
#       data_mod$data2$kriging_rss$RSID <- D2$snp
#       data_mod$data2$kriging_rss$outlier <- ifelse(data_mod$data2$kriging_rss$logLR  > 2 &
#                                                      abs(data_mod$data2$kriging_rss$z) > 2, TRUE, FALSE)
#       data_mod$source <- "Dataset 2"
#       data_mod$data   <- calc_credible_set(rsids          = names(data_mod$data2$pip),
#                                            snp.pp         = data_mod$data2$pip,
#                                            dat_join_to    = data.table::copy(source_2_mod$data),
#                                            credible_sets  = data_mod$data2$sets$cs,
#                                            credible_set_p = NULL)
#
#
#       # DATASET 1&2 COLOC
#     } else if(input$method=="Coloc" && input$assumption == "Single" && !is.null(source_1_mod$data) && !is.null(source_2_mod$data)) {
#
#       # harmonise the datasets
#       harm <- genepi.utils::harmonise(gwas1 = source_1_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
#                                       gwas2 = source_2_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
#                                       gwas1_trait = "d1",
#                                       gwas2_trait = "d2",
#                                       merge = c("SNP"="SNP"))
#
#       # make and check the dataset 1
#       D1 <- make_coloc_dataset(dat  = harm[, list(RSID=RSID_d1,BP=BP_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
#                                type = input$source_1_type,
#                                sdY  = input$sd_y1,
#                                ld   = NULL)
#
#       # make and check the dataset 2
#       D2 <- make_coloc_dataset(dat  = harm[, list(RSID=RSID_d2,BP=BP_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
#                                type = input$source_2_type,
#                                sdY  = input$sd_y2,
#                                ld   = NULL)
#
#       # run the coloc abf function
#       res <- coloc::coloc.abf(D1, D2,
#                               p1  = input$coloc_p1,
#                               p2  = input$coloc_p2,
#                               p12 = input$coloc_p12)
#
#       # assign the results
#       data_mod$data2 <- res
#       data_mod$source <- "Coloc"
#
#       # populate data_mod$data with whichever dataset we plan on feeding downstream
#       if(input$downstream_dataset == "Dataset 1") {
#
#         dat_join_to = source_1_mod$data
#
#       } else if(input$downstream_dataset == "Dataset 2") {
#
#         dat_join_to = source_2_mod$data
#
#       }
#       data_mod$data <- calc_credible_set(rsids          = data_mod$data2$results$snp,
#                                          snp.pp         = data_mod$data2$results$SNP.PP.H4,
#                                          dat_join_to    = dat_join_to,
#                                          credible_set_p = input$credible_set_p)
#
#       # DATASET 1&2 Coloc SuSiE
#     } else if(input$method=="Coloc" && input$assumption == "SuSiE" && !is.null(source_1_mod$data) && !is.null(source_2_mod$data)) {
#
#       # harmonise the datasets
#       harm <- genepi.utils::harmonise(gwas1 = source_1_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
#                                       gwas2 = source_2_mod$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,N)],
#                                       gwas1_trait = "d1",
#                                       gwas2_trait = "d2",
#                                       merge = c("SNP"="SNP"))
#
#       # get the reference file
#       plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
#                                    chrom    = app$modules$gene$chr,
#                                    from     = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
#                                    to       = app$modules$gene$end + app$modules$gene$flanks_kb*1000)
#
#       # create the LD matrix for the region variants
#       dat_ld_obj1 <- genepi.utils::ld_matrix(variants     = harm[, list(RSID=RSID_d1,BP=BP_d1,EA=EA_d1,OA=OA_d1,BETA=BETA_d1,SE=SE_d1,EAF=EAF_d1,P=P_d1,N=N_d1)],
#                                              with_alleles = TRUE,
#                                              plink2       = get_plink2_exe(),
#                                              plink_ref    = plink_ref)
#       dat_ld_obj2 <- genepi.utils::ld_matrix(variants     = harm[, list(RSID=RSID_d2,BP=BP_d2,EA=EA_d2,OA=OA_d2,BETA=BETA_d2,SE=SE_d2,EAF=EAF_d2,P=P_d2,N=N_d2)],
#                                              with_alleles = TRUE,
#                                              plink2       = get_plink2_exe(),
#                                              plink_ref    = plink_ref)
#
#       # make and check the dataset 1
#       D1 <- make_coloc_dataset(dat  = dat_ld_obj1[["dat"]],
#                                type = input$source_1_type,
#                                sdY  = input$sd_y1,
#                                ld   = dat_ld_obj1[["ld_mat"]])
#
#       # make and check the dataset 2
#       D2 <- make_coloc_dataset(dat  = dat_ld_obj2[["dat"]],
#                                type = input$source_2_type,
#                                sdY  = input$sd_y2,
#                                ld   = dat_ld_obj2[["ld_mat"]])
#
#
#       # run SuSiE
#       D1_susie <- coloc::runsusie(D1, coverage=input$credible_set_p)
#       D2_susie <- coloc::runsusie(D2, coverage=input$credible_set_p)
#
#       # run the coloc abf function
#       res <- coloc::coloc.susie(D1_susie, D2_susie)
#
#
#       # Browse[1]>  print(res$summary)
#       # nsnps            hit1            hit2     PP.H0.abf     PP.H1.abf     PP.H2.abf    PP.H3.abf     PP.H4.abf  idx1  idx2
#       # <int>          <char>          <char>         <num>         <num>         <num>        <num>         <num> <int> <int>
#       #   1:   390  rs10305435_C_T  rs10305435_C_T 2.190624e-175  2.093143e-89  2.093143e-89 0.000000e+00  1.000000e+00     1     1
#       # 2:   390  rs10305445_C_T  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  1.628576e-88     2     1
#       # 3:   390   rs2300613_G_A  rs10305435_C_T 5.057475e-161  1.046572e-86  4.832421e-75 1.000000e+00  9.458692e-46     3     1
#       # 4:   390   rs2300614_A_G  rs10305435_C_T 6.092174e-263  1.046572e-86 5.821077e-177 1.000000e+00  8.832468e-88     4     1
#       # 5:   390    rs877446_A_G  rs10305435_C_T 3.379615e-298  1.046572e-86 3.229225e-212 1.000000e+00  1.935952e-85     6     1
#       # 6:   390  rs10305442_A_G  rs10305435_C_T 2.369280e-277  1.046572e-86 2.263849e-191 1.000000e+00  1.214910e-88     7     1
#       # 7:   390 rs115493670_C_T  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  4.568936e-89     8     1
#       # 8:   390   rs9283907_A_G  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  1.434933e-88     9     1
#       # 9:   390   rs9283906_T_C  rs10305435_C_T  0.000000e+00  1.046572e-86  0.000000e+00 1.000000e+00  4.752047e-89    10     1
#       # 10:   390   rs2235868_A_C  rs10305435_C_T 3.821291e-149  1.046572e-86  3.651246e-63 1.000000e+00  3.825684e-25     5     1
#       # 11:   390  rs10305435_C_T  rs10305445_C_T  0.000000e+00  0.000000e+00  1.046572e-86 1.000000e+00  1.628576e-88     1     2
#
#       #
#       #
#       # # assign the results
#       # data_mod$data2 <- res
#       # data_mod$source <- "Coloc"
#       #
#       # # populate data_mod$data with whichever dataset we plan on feeding downstream
#       # if(input$downstream_dataset == "Dataset 1") {
#       #
#       #   dat_join_to = source_1_mod$data
#       #
#       # } else if(input$downstream_dataset == "Dataset 2") {
#       #
#       #   dat_join_to = source_2_mod$data
#       #
#       # }
#       # data_mod$data <- calc_credible_set(rsids          = data_mod$data2$results$snp,
#       #                                    snp.pp         = data_mod$data2$results$SNP.PP.H4,
#       #                                    dat_join_to    = dat_join_to,
#       #                                    credible_set_p = input$credible_set_p)
#       #
#
#     }
#
#     # button text back to black
#     shinyjs::runjs(paste0('document.getElementById("', ns("run"), '").style.color = "black";'))
#
#   },
#   error=function(e) {
#     showNotification(paste0("coloc::finemap.abf() failed - ", e), type="error", duration=10)
#     return(NULL)
#   })
#
#
# })
#
#
# #==========================================
# # Observe controls - make run button red if changed
# #==========================================
# observeEvent(list(input$downstream_dataset,
#                   input$source_1,
#                   input$source_2,
#                   input$source_1_type,
#                   input$source_2_type,
#                   input$sd_y1,
#                   input$sd_y2,
#                   input$coloc_p1,
#                   input$coloc_p12,
#                   input$coloc_p2,
#                   input$coloc_h4,
#                   input$ld_reference,
#                   input$assumption,
#                   input$method,
#                   input$credible_set_p
# ), {
#
#   shinyjs::runjs(paste0('document.getElementById("', ns("run"), '").style.color = "red";'))
#
# }, ignoreInit=TRUE)
#
#
# #==========================================
# # Observe credible set P-value threshold
# #==========================================
# session$userData[[ns("credible_set_p")]] <- observeEvent(input$credible_set_p, {
#
#   if(is.null(data_mod$data)) return(NULL)
#   req("cumsum_pp_h4" %in% names(data_mod$data))
#
#   # determine credible set - flag as coloc
#   data_mod$data[, coloc := ifelse(cumsum_pp_h4 > input$credible_set_p, TRUE, FALSE)]
#
# })
#
#
# #==========================================
# # Observe additions / deletions of modules
# #==========================================
# session$userData[[ns("app-modules")]] <- observeEvent(app$modules, {
#   updateSelectInput(session, inputId="source_1", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
#   updateSelectInput(session, inputId="source_2", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
# })
#
#
# #==========================================
# # Observe input selection for choices init as observing just app$modules doesn't work at module init...
# #==========================================
# session$userData[[ns("sources")]] <- observeEvent(list(input$source_1,input$source_2), {
#   if(is.null(input$source_1) || input$source_1=="") {
#     updateSelectInput(session, inputId="source_1", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
#   }
#   if(is.null(input$source_2) || input$source_2=="") {
#     updateSelectInput(session, inputId="source_2", choices=names(app$modules)[!names(app$modules)  %in% c("gene",id)])
#   }
# })
#
#
# #==========================================
# # Observe method input (SuSiE / single causal variant assumption)
# #==========================================
# session$userData[[ns("assumption")]] <- observeEvent(input$assumption, {
#   if(input$assumption == "SuSiE") {
#     shinyjs::enable("ld_reference")
#   } else if(input$assumption == "Single") {
#     shinyjs::disable("ld_reference")
#   }
# })



