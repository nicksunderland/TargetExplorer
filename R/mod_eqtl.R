#' eQTL UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_eqtl_ui <- function(id){
  ns <- NS(id)
  # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
  div(
    id = id,
    sidebarLayout(position = "right",
                  sidebarPanel(p(strong(paste0("Controls [",id,"]"))),
                               width = 3,
                               hr(),
                               mod_data_ui(id=paste0(id,"-data")),
                               hr(),
                               mod_clump_ui(id=ns("clump")),
                               hr(),
                               fluidRow(
                                 column(6,
                                        selectInput(inputId = ns("source_2"),
                                                    label   = "Upstream data",
                                                    choices = c(""))
                                 ),
                                 column(6,
                                        selectInput(inputId = ns("source_2_filter"),
                                                    label   = "Upstream filter",
                                                    choices = c("All"),
                                                    selected = "All")
                                 )
                               ),
                               fluidRow(
                                 column(6,
                                        selectizeInput(inputId = ns("datasets"),
                                                       label   = "Datasets",
                                                       choices = c("All"),
                                                       selected= "All",
                                                       multiple= TRUE),
                                        prettyRadioButtons(inputId = ns("beta_dir"),
                                                           label   = "BETA dir",
                                                           choices = c("Any","Concordant","Disconcordant"),
                                                           selected = "Any",
                                                           inline   = FALSE)
                                 ),
                                 column(6,
                                        selectInput(inputId = ns("dataset_combine_method"),
                                                    label   = "Combine by",
                                                    choices = c("Don't combine", "Maximum abs(BETA)","Lowest P-value","Largest sample size"),
                                                    selected = "Largest sample size"),
                                        actionButton(inputId = ns("apply"),
                                                     label   = "Apply"),
                                        actionButton(inputId = ns("reset"),
                                                     width   = "40px",
                                                     label   = "",
                                                     icon    = icon("rotate-left"))
                                 )
                               ),

                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            fluidRow(
                              column(1, p(strong("eQTL:"))),
                              column(4, selectInput(inputId = ns("view"),
                                                    label   = NULL,
                                                    choices = c("Tissue - variation", "Tissue - expression", "Tissue - clusters"),
                                                    selected = "Tissue overview")),
                            ),
                            fluidRow(
                              column(6,
                                     plotOutput(outputId = ns("tissue_plot"),
                                                height   = "500px",
                                                brush    = ns("mr_plot_brush")),
                                     tableOutput(outputId = ns("mr_result"))
                              ),
                              column(6,
                                     plotOutput(outputId = ns("variants_plot"),
                                                height   = "500px",
                                                brush    = ns("locus_plot_brush"))
                              ),
                              tableOutput(outputId = ns("locus_plot_table"))
                            ),
                  ) # main panel end
    ), # sidebar layout end
    hr()
  ) # div end
}


#' eQTL Server Functions
#' @import ggplot2 ggrepel ggfortify
#' @importFrom stats dnorm sd
#' @noRd
mod_eqtl_server <- function(id, app){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # R CMD checks
    GENE_NAME <- BP_END <- BP_START <- GENE_NAME <- nlog10P <- BETA_source <- BETA_eqtl <- EAF <- EA <- OA <- CHR <- BP <- N <- P <- eqtl <- NULL
    DATASET <- TPM <- mean_tpm <- QUANT_METHOD <- STUDY_TISSUE_METHOD <- TISSUE <- STUDY <- STUDY_TISSUE <- RSID <- RSID_beta <- CHR <- BP <- N <- P <- eqtl <- NULL
    sd_beta <- BETA <- mean_beta <- beta_dnorm <- NULL

    #==========================================
    # Data module server for the eQTL module
    #==========================================
    data_mod  <- mod_data_server(id="data", gene_module=app$modules$gene)
    clump_mod <- mod_clump_server(id="clump", gene_module=app$modules$gene, data_module=data_mod)


    #==========================================
    # Observe additions / deletions of upstream source modules
    #==========================================
    observeEvent(app$modules, {
      updateSelectInput(session, inputId="source_2", choices=c("", names(app$modules)[!names(app$modules) %in% c("gene",id)]))
    })


    #==========================================
    # Observe input selection for choices init as observing just app$modules doesn't work at module init...
    #==========================================
    observeEvent(input$source_2, {
      if(is.null(input$source_2) || input$source_2=="") {
        updateSelectInput(session, inputId="source_2", choices=c("", names(app$modules)[!names(app$modules) %in% c("gene",id)]))
      }
      if(!is.null(app$modules[[input$source_2]]$data)) {
        source_cols <- names(app$modules[[input$source_2]]$data)
        possible_filters <- source_cols[source_cols %in% c('index','coloc')]
        updateSelectInput(session, inputId="source_2_filter", choices=c("All",possible_filters))
      }
    })


    #==========================================
    # Process input data when it changes / is loaded; update GUI also
    #==========================================
    observeEvent(data_mod$data, {
      req(data_mod$data)

      # mean beta per quantification method
      data_mod$data[, mean_beta := mean(BETA, na.rm=TRUE), by="QUANT_METHOD"]
      data_mod$data[, sd_beta   := stats::sd(BETA, na.rm=TRUE),   by="QUANT_METHOD"]

      # beta dnorm and ordered RSID factor (nicer plotting)
      data_mod$data[, beta_dnorm := stats::dnorm(BETA, mean=mean_beta, sd=sd_beta)]
      data.table::setorder(data_mod$data, beta_dnorm)
      data_mod$data[, RSID_beta := factor(RSID, levels=unique(data_mod$data$RSID))]

      # Tissue - study factor
      data_mod$data[, STUDY_TISSUE := factor(paste(STUDY,"-",TISSUE))]
      data_mod$data[, STUDY_TISSUE_METHOD := factor(paste(STUDY,"-",TISSUE,"-",QUANT_METHOD))]

      # mean expression
      data_mod$data[, mean_tpm  := mean(TPM, na.rm=TRUE),  by="STUDY_TISSUE"]

      # transcripts per million dnorm and ordered RSID factor (nicer plotting)
      data.table::setorder(data_mod$data, mean_tpm)
      data_mod$data[, TISSUE := factor(TISSUE, levels=unique(data_mod$data$TISSUE))]

      # Calculate PCA / clustering datatable
      beta_matrix <- data_mod$data[QUANT_METHOD=="ge" & grepl("^rs[0-9]+$",RSID), list(RSID, DATASET=paste(STUDY,"-",TISSUE), BETA)]
      beta_matrix <- beta_matrix[, .SD[which.max(BETA)], by=c("RSID","DATASET")] # a few duplicates in a dataset, remove
      beta_matrix <- tidyr::complete(beta_matrix, RSID, DATASET) |> data.table::as.data.table()
      beta_matrix <- data.table::dcast(beta_matrix, RSID ~ DATASET, value.var="BETA")
      beta_matrix <- stats::na.omit(beta_matrix)
      data_mod$pca_data <- data.table::as.data.table( t(beta_matrix[,-1]) )
      data.table::setnames(data_mod$pca_data, names(data_mod$pca_data), beta_matrix$RSID)
      data_mod$pca_data[, DATASET := names(beta_matrix)[-1]]
      data.table::setcolorder(data_mod$pca_data, "DATASET")

      # update filtering options
      current_datasets <- input$datasets
      updateSelectInput(session, inputId="datasets", choices=c("All", unique(levels(data_mod$data$STUDY_TISSUE_METHOD))), selected=current_datasets)

    })


    #==========================================
    # eQTL tissue / database analysis plots
    #==========================================
    output$tissue_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(data_mod$data), 'No data imported, click the import button')
      )

      if(input$view=="Tissue - variation") {

        p <- ggplot2::ggplot(data = data_mod$data,
                             mapping = ggplot2::aes(x = RSID_beta, y = STUDY, fill=beta_dnorm)) +
          ggplot2::geom_tile() +
          viridis::scale_fill_viridis(option="inferno", direction=-1) +
          ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
          ggplot2::labs(x = "RSID", y = "Dataset (study)", fill="dnorm(\u03B2)") +
          facet_grid(rows = vars(TISSUE), cols = vars(QUANT_METHOD))

      } else if(input$view=="Tissue - expression") {

        p <- ggplot2::ggplot(data = data_mod$data,
                             mapping = ggplot2::aes(x = TISSUE, y = STUDY, fill=mean_tpm)) +
          ggplot2::geom_tile() +
          viridis::scale_fill_viridis(option="inferno", direction=1) +
          ggplot2::labs(x = "Tissue", y = "Study", fill="TPM")

      } else if(input$view=="Tissue - clusters") {

        # check data
        validate(
          need(nrow(data_mod$pca_data)>1, 'No PCA analysis avilable - only one tissue/dataset')
        )

        # run PCA
        pca <- stats::prcomp(data_mod$pca_data[,-1], center=TRUE, scale=TRUE)

        # plot
        p <- autoplot(pca, data=data_mod$pca_data, color='DATASET') +
          geom_point(aes(color=DATASET), size=5) +
          theme_classic() +
          labs(color="Dataset")
      }

      return(p)

    })


    #==========================================
    # reset filtering of upstream variants
    #==========================================
    observeEvent(input$reset, {

      if("eqtl" %in% names(data_mod$data)) {
        data_mod$data[, eqtl := NULL]
        tmp <- data.table::copy(data_mod$data)
        data_mod$data <- NULL
        data_mod$data <- tmp
      }

    })

    #==========================================
    # apply filtering of upstream variants
    #==========================================
    observeEvent(input$apply, {

      # require the inputs (i.e. not NULL)
      req(data_mod$data, input$datasets, input$beta_dir, input$dataset_combine_method, input$source_2, input$source_2_filter)

      # which eQTL dataset to use?
      if("All" %in% input$datasets) {

        data_mod$data[, eqtl := TRUE]

      } else {

        data_mod$data[, eqtl := ifelse(STUDY_TISSUE_METHOD %in% input$datasets,
                                       rep(TRUE, nrow(data_mod$data)),
                                       rep(FALSE, nrow(data_mod$data)))]

        if(input$dataset_combine_method == "Maximum abs(BETA)") {

          data_mod$data[eqtl==TRUE, eqtl := ifelse(abs(BETA) == max(abs(BETA), na.rm=T) & !duplicated(abs(BETA) == max(abs(BETA), na.rm=T)), # ensure only one TRUE per RSID (BETA==max might give >1 value)
                                                   rep(TRUE,.N), rep(FALSE,.N)), by=RSID]

        } else if(input$dataset_combine_method == "Lowest P-value") {

          data_mod$data[eqtl==TRUE, eqtl := ifelse(P == min(P, na.rm=T) & !duplicated(P == min(P, na.rm=T)), # ensure only one TRUE per RSID (P==min(P) might give >1 value)
                                                   rep(TRUE,.N), rep(FALSE,.N)), by=RSID]

        } else if(input$dataset_combine_method == "Largest sample size") {

          data_mod$data[eqtl==TRUE, eqtl := ifelse(N == max(N, na.rm=T) & !duplicated(N == max(N, na.rm=T)), # ensure only one TRUE per RSID (N==max(N) might give >1 value)
                                                   rep(TRUE,.N), rep(FALSE,.N)), by=RSID]

        }

      }

      # which variants are also in the upstream source data
      if(input$source_2 != "" && !is.null(app$modules[[input$source_2]]$data)) {

        if(input$source_2_filter != "All") {

          # all RSIDs where the filter condition in the upstream data is TRUE
          upstream_variants <- app$modules[[input$source_2]]$data$RSID[ which(app$modules[[input$source_2]]$data[[input$source_2_filter]]==TRUE) ]

        } else {

          # all RSIDs in the upstream data
          upstream_variants <- app$modules[[input$source_2]]$data$RSID

        }

        # find the corresponding variants in the eQTL data; set false any not in the upstream, else carry the eqtl flag forward
        data_mod$data[, eqtl := ifelse(!RSID %in% upstream_variants, FALSE, eqtl)]

      }

      # beta directionality - requires harmonisation
      if(input$source_2 != "" && !is.null(app$modules[[input$source_2]]$data) && input$beta_dir != "Any") {

        # relevant eqtl data and source data so far
        eqtl_data   <- data_mod$data[eqtl==TRUE, list(SNP=RSID,CHR,BP,EA,OA,EAF,BETA,P)]
        source_data <- app$modules[[input$source_2]]$data[RSID %in% eqtl_data$SNP, list(SNP=RSID, CHR,BP,EA,OA,EAF,BETA,P)]

        # harmonise the data
        harm_dat <- genepi.utils::harmonise(gwas1 = eqtl_data,
                                            gwas2 = source_data,
                                            gwas1_trait = "eqtl",
                                            gwas2_trait = "source",
                                            merge = c("SNP"="SNP"))

        # which direction?
        if(input$beta_dir == "Concordant") {

          harm_dat <- harm_dat[sign(BETA_eqtl) == sign(BETA_source), ]

        } else if(input$beta_dir == "Disconcordant") {

          harm_dat <- harm_dat[sign(BETA_eqtl) != sign(BETA_source), ]

        }

        # flag which ones are left, else carry the eqtl flag forward (allow palindromic SNPs to carry forward)
        data_mod$data[, eqtl := ifelse(!RSID %in% harm_dat$SNP_eqtl[harm_dat$keep | harm_dat$palindromic] , FALSE, eqtl)]

      }

      # force memory location change (reactivity doesnt work with inplace data.table modifications)
      tmp <- data.table::copy(data_mod$data)
      data_mod$data <- NULL
      data_mod$data <- tmp
    })

    #==========================================
    # eQTL chosen / filtered variants plot
    #==========================================
    output$variants_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(data_mod$data), 'No data imported, click the import button')
      )

      print("plotting")

      p <- ggplot2::ggplot(data = data_mod$data,
                           mapping = ggplot2::aes(x = BP, y = nlog10P, color=STUDY_TISSUE, shape=QUANT_METHOD)) +
        ggplot2::geom_point() +
        annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        viridis::scale_color_viridis(option="viridis", discrete = TRUE) +
        lims(x = c(app$modules$gene$start - app$modules$gene$flanks_kb*1000, app$modules$gene$end + app$modules$gene$flanks_kb*1000),
             y = c(-2.5, max(8.0, max(data_mod$data$nlog10P)))) +
        labs(x     = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y     = expression(paste("-log"[10], plain(P))),
             color = "Study - Tissue",
             shape = "Quant. method")

      # if gene data then add
      if(!is.null(data_mod$data2)) {
        p <- p +
          geom_hline(yintercept = 0) +
          geom_rect(data = data_mod$data2[data_mod$data2$STRAND=="+",  ],
                    inherit.aes=FALSE,
                    mapping = aes(xmin=BP_START, xmax=BP_END, ymin=-1.75, ymax=-1.25), fill="grey", alpha=0.3) +
          geom_rect(data = data_mod$data2[data_mod$data2$STRAND=="-",  ],
                    inherit.aes=FALSE,
                    mapping = aes(xmin=BP_START, xmax=BP_END, ymin=-2.5, ymax=-2), fill="grey", alpha=0.3) +
          geom_text_repel(data = data_mod$data2[data_mod$data2$STRAND=="-",  ],
                          mapping = aes(label = GENE_NAME, x=(BP_END-BP_START)/2 + BP_START, y =-2.25),
                          inherit.aes=FALSE, min.segment.length = 0.25) +
          geom_text_repel(data = data_mod$data2[data_mod$data2$STRAND=="+",  ],
                          mapping = aes(label = GENE_NAME, x=(BP_END-BP_START)/2 + BP_START, y =-1.5),
                          inherit.aes=FALSE, min.segment.length = 0.25)
      }


      # if we have computed the eQTL flag then plot
      if("eqtl" %in% names(data_mod$data)) {

        # draw large triangles, lines and labels if less than 50 to draw; otherwise small triangles
        if(sum(data_mod$data$eqtl, na.rm=TRUE) < 50) {

          p <- p +
            geom_vline(data = data_mod$data[eqtl==TRUE, ], mapping=aes(xintercept=BP), linetype="dotted", color= "darkred") +
            geom_point(data = data_mod$data[eqtl==TRUE, ], mapping=aes(x=BP, y=nlog10P), inherit.aes=FALSE, size=3, fill="red",color="red",shape=24) +
            geom_label_repel(data = data_mod$data[eqtl==TRUE, ], mapping=aes(label=RSID, x=BP, y=nlog10P), max.overlaps=Inf, inherit.aes=FALSE) +
            labs(color = "Upstream hit", fill = "Upstream hit")

        } else {

          p <- p +
            geom_point(data = data_mod$data[eqtl==TRUE, ], mapping=aes(x=BP, y=nlog10P), inherit.aes=FALSE, size=1, fill="red",color="red",shape=24) +
            labs(color = "Upstream hit", fill = "Upstream hit")

        }

      }

      return(p)

    })


    #==========================================
    # Point select table - for brushed point
    #==========================================
    output$locus_plot_table <- renderTable({
      tryCatch({
        bp <- brushedPoints(data_mod$data, input$locus_plot_brush)
        if(nrow(bp)==0 && !is.null(data_mod$data2)) {
          bp <- data_mod$data2[(input$locus_plot_brush$xmin < data_mod$data2$BP_START) & (input$locus_plot_brush$xmax > data_mod$data2$BP_END) &
                                 ((input$locus_plot_brush$ymin < -1.75 & input$locus_plot_brush$ymax > -1.25 & data_mod$data2$STRAND=="+") |
                                    (input$locus_plot_brush$ymin < -2.50 & input$locus_plot_brush$ymax > -2.00 & data_mod$data2$STRAND=="-")), ]
        }
        if(nrow(bp)>0) return(bp) else return(NULL)
      },
      error = function(e) {
        return(NULL)
      })
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
