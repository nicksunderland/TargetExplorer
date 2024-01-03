#' GWAS UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#'
mod_gwas_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_gwas_ui(ns=",ns("foo"),")"))
  # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
  div(
    id = id,
    sidebarLayout(position = "right",
                  sidebarPanel(width = 3,
                               fluidRow(column(10, p(strong(paste0("Controls [",id,"]")))),
                                        column(2,  mod_remove_ui(id=ns("remove")))),
                               hr(),
                               mod_data_ui(id=ns("data")),
                               hr(),
                               mod_grouping_ui(id=ns("grouping"))
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            fluidRow(
                              column(6, p(strong("Locus plot:"))),
                              column(1, p(strong("Sensitivity:"), style="text-align:right; margin-top: 6px;")),
                              column(2, selectInput(inputId = ns("sensitivity_plot"), label=NULL, choices = c("Off"), selected = "Off")),
                              column(1, p(strong("Source:"), style="text-align:right; margin-top: 6px;")),
                              column(2, mod_source_select_ui(id=ns("sensitivity_source")))
                            ),
                            # Locus zoom plot
                            uiOutput(outputId = ns("plot_area")),
                            tableOutput(outputId = ns("locus_plot_table")),
                  ) # main panel end
    ), # sidebar layout end
    hr()
  ) # div end
}


#' gwas Server Functions
#' @import ggplot2 ggrepel
#' @noRd
mod_gwas_server <- function(id, app){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_gwas_server(ns=",ns("foo"),")"))

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- group <- log10P <- SNP <- index <- nlog10P <- NULL


    #==========================================
    # Data, clump, and remove (sub)module servers for the GWAS module
    #==========================================
    data_mod     <- mod_data_server(id="data", gene_module=app$modules$gene)
    grouping_mod <- mod_grouping_server(id="grouping", gene_module=app$modules$gene, data_module=data_mod)
    remove_mod   <- mod_remove_server(id="remove", app=app, parent_id=id, parent_inputs=input)
    sensitivity_source_mod <- mod_source_select_server(id="sensitivity_source", app=app, source_type=c("GWAS","eQTL","Coloc"), label=NULL)


    #==========================================
    # Observe the grouping input
    #==========================================
    observeEvent(list(grouping_mod$grouping, data_mod$source), {

      req(grouping_mod$grouping)
      req(data_mod$source)

      source_specific_plots <- switch(data_mod$source,
                                      "Local"              = c(),
                                      "IEU Open GWAS"      = c(),
                                      "EBI GWAS Catalogue" = c(),
                                      "EBI eQTL Catalogue" = c("Tissue - variation","Tissue - expression","Tissue - clusters"))

      sensitivity_plot_choices <- switch(grouping_mod$grouping,
                                         "Off"          = c("Off", source_specific_plots),
                                         "plink::clump" = c("Off", "LD structure", "Clump exp-out corr", source_specific_plots),
                                         "r-coloc"      = c("Off", "LD structure", "Kriging plot", source_specific_plots),
                                         "r-susieR"     = c("Off", "LD structure", "Kriging plot", source_specific_plots))

      updateSelectInput(inputId = "sensitivity_plot", choices = sensitivity_plot_choices)

    })


    #==========================================
    # Observe the sensitivity input
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
    # Manhattan / LocusZoom plot
    #==========================================
    output$locus_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(data_mod$data), 'No data imported, click the import button')
      )

      # create the locus plot
      p <- ggplot(data    = data_mod$data,
                  mapping = aes(x=BP, y=nlog10P)) +
        geom_hline(yintercept=7.3, linetype="dotted", color="lightgrey") +
        annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        lims(x = c(app$modules$gene$start - app$modules$gene$flanks_kb*1000, app$modules$gene$end + app$modules$gene$flanks_kb*1000),
             y = c(-2.5, max(8.0, max(data_mod$data$nlog10P)))) +
        labs(x        = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y        = expression(paste("-log"[10], plain(P))),
             subtitle = data_mod$source_id)

      # if alternative color/group data (and no grouping column from the group module), plot that color
      possible_base_grouping <- list(eqtl = c("STUDY","TISSUE","QUANT_METHOD"))
      if(!all(c("index","group") %in% names(data_mod$data)) &&
         any(sapply(possible_base_grouping, function(cols) all(cols %in% names(data_mod$data))))) {

        # must have been imported from eQTL catalogue
        if(all(possible_base_grouping[["eqtl"]] %in% names(data_mod$data))) {

          p <- p +
            ggplot2::geom_point(aes(color=TISSUE, shape=QUANT_METHOD)) +
            labs(color  = "Tissue",
                 shape = "Quantification")

        }

      # else just color grey
      } else {

        p <- p + geom_point(color="lightgray")

      }

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
                          inherit.aes=FALSE,
                          direction = "x", min.segment.length = 0.25) +
          geom_text_repel(data = data_mod$data2[data_mod$data2$STRAND=="+",  ],
                          mapping = aes(label = GENE_NAME, x=(BP_END-BP_START)/2 + BP_START, y =-1.5),
                          inherit.aes=FALSE,
                          direction = "x", min.segment.length = 0.25)
      }

      # if there is grouped (clumped or finemapped) data, plot
      if(all(c("index","group") %in% colnames(data_mod$data))) {
        p <- p +
          geom_point(data = data_mod$data[!is.na(data_mod$data$group), ],            mapping = aes(x=BP, y=nlog10P, color=group, fill=group), shape=23) +
          geom_vline(data = data_mod$data[which(data_mod$data$index==TRUE), ],       mapping = aes(xintercept=BP), linetype="dotted", color="darkred") +
          geom_point(data = data_mod$data[which(data_mod$data$index==TRUE), ],       mapping = aes(x=BP, y=nlog10P), size=3, fill="red", color="red", shape=24) +
          geom_label(data = data_mod$data[which(data_mod$data$index==TRUE), ],       mapping = aes(label=group, x=BP, y=-0.5)) +
          geom_label_repel(data = data_mod$data[which(data_mod$data$index==TRUE), ], mapping = aes(label=RSID,  x=BP, y=nlog10P)) +
          labs(color = "Group", fill = "Group")


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
    # Sensitivity plots
    #==========================================
    output$sensitivity_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(input$sensitivity_plot) && input$sensitivity_plot != "Off", 'No sensitivity plots selected'),
        need(!is.null(data_mod$data), 'No GWAS data loaded')
      )


      #--------------------------------------------
      # LD structure plots
      if(input$sensitivity_plot == "LD structure") {

        validate(
          need(!is.null(data_mod$ld_matrix_r2) | !is.null(data_mod$ld_matrix_r), 'No LD data - run clumping'),
        )

        # R2 data
        if(is.null(data_mod$ld_matrix_r2)) {
          p_r2 <- NULL
        } else {
          upper_tri <- data_mod$ld_matrix_r2
          upper_tri[lower.tri(upper_tri)] <- NA
          upper_tri <- reshape2::melt(upper_tri, na.rm=TRUE)
          p_r2 <- ggplot(data = upper_tri,
                         mapping = aes(Var1, Var2, fill=value)) +
            geom_tile() +
            scale_fill_gradient(low="lightyellow", high="red", limits=c(0,1)) +
            labs(fill = "r2") +
            theme(axis.title.x = element_blank(),
                  axis.text.x  = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y  = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = c(.95, .05),
                  legend.justification = c("right", "bottom"),
                  legend.box.just = "right",
                  legend.margin = margin(6, 6, 6, 6))
        }

        # R corr data
        if(is.null(data_mod$ld_matrix_r)) {
          p_r <- NULL
        } else {
          upper_tri <- data_mod$ld_matrix_r
          upper_tri[lower.tri(upper_tri)] <- NA
          upper_tri <- reshape2::melt(upper_tri, na.rm=TRUE)
          p_r  <- ggplot(data = upper_tri,
                         mapping = aes(Var1, Var2, fill=value)) +
            geom_tile() +
            scale_fill_gradient2(low="darkblue", mid = "white", high="red", midpoint = 0, limits=c(-1,1)) +
            labs(fill = "r") +
            theme(axis.title.x = element_blank(),
                  axis.text.x  = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y  = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = c(.95, .05),
                  legend.justification = c("right", "bottom"),
                  legend.box.just = "right",
                  legend.margin = margin(6, 6, 6, 6))
        }

        # combine the heat maps
        p <- ggpubr::ggarrange(plotlist = list(p_r2, p_r), ncol = 1)


     #--------------------------------------------
     # Clumping exposure-outcome correlation plots
     } else if(input$sensitivity_plot == "Clump exp-out corr") {

       validate(
         need("group" %in% names(data_mod$data), 'No clumping data found - run clumping'),
         need(!is.null(sensitivity_source_mod$data), 'No source data - select a source'),
       )

       # harmonise the datasets
       data_mod$data[, SNP:= RSID]
       sensitivity_source_mod$data[, SNP:= RSID]
       harm <- genepi.utils::harmonise(data_mod$data, sensitivity_source_mod$data, gwas1_trait="exposure", gwas2_trait="outcome", merge = c("RSID"="RSID"))
       harm <- harm[keep==TRUE, ]
       if("group_exposure" %in% names(harm)) {
         data.table::setnames(harm, "group_exposure", "group")
       }

       # plot
       p <- ggplot(data    = harm[!is.na(group), ],
                   mapping = aes(x=BETA_exposure, y=BETA_outcome, color=group)) + #alpha=P_exposure
         geom_point() +
         # scale_alpha_discrete(name="Exposure P-value") +
         geom_smooth(method="lm", mapping = aes(weight = (1/SE_exposure)*(1/SE_outcome)), color="red") +
         geom_point(data = harm[which(harm$index==TRUE), ], mapping = aes(x=BETA_exposure, y=BETA_outcome), size=3, fill="red", color="red", shape=24) +
         # theme_classic() +
         labs(x = '\u03B2 - exposure',
              y = '\u03B2 - outcome',
              color = "Group") +
         facet_wrap(~group, nrow=4, scales = "free")


     #--------------------------------------------
     # Kriging plot for LD reference vs study fit
     } else if(input$sensitivity_plot == "Kriging plot") {

       validate(
         need(!is.null(data_mod$kriging_rss), 'No Kriging plot data - run finemapping'),
       )

       p <- ggplot(data    = data_mod$kriging_rss,
                   mapping = aes(x=condmean, y=z)) +
         geom_abline(intercept = 0, slope = 1) +
         geom_point(mapping = aes(color=outlier)) +
         geom_label_repel(data = data_mod$kriging_rss[data_mod$kriging_rss$outlier==TRUE, ],
                          mapping = aes(label=RSID), max.overlaps = Inf) +
         scale_color_manual(values = c("TRUE"="red", "FALSE"="darkgrey"), labels=c("TRUE"="Outlier","FALSE"="Within tolerance"), drop=FALSE) +
         theme_classic() +
         labs(y = "Observed z scores", x = "Expected value") +
         theme(legend.position = "top",
               legend.title = element_blank())


     #--------------------------------------------
     # Tissue variation (betas of different tissues) plot - for eQTL data
     } else if(input$sensitivity_plot == "Tissue - variation") {

       # order by beta for nicer plotting
       data.table::setorder(data_mod$data, BETA)
       data_mod$data[, RSID_beta := factor(RSID, levels=unique(data_mod$data$RSID))]

       # plot
       p <- ggplot(data    = data_mod$data,
                   mapping = aes(x = RSID_beta, y = STUDY, fill=BETA)) +
         geom_tile() +
         scale_fill_gradient2(low="red",mid="black",high="green", limits=c(-0.5,0.5), oob=scales::squish) +
         theme(axis.text.x = element_blank()) +
         labs(x = "RSID", y = "Dataset (study)", fill="\u03B2") +
         facet_grid(rows = vars(TISSUE), cols = vars(QUANT_METHOD))


     #--------------------------------------------
     # Tissue expression (TPM of different tissues) plot - for eQTL data
     } else if(input$sensitivity_plot == "Tissue - expression") {

       # order by transcripts per million (nicer plotting)
       data.table::setorder(data_mod$data, TPM)
       data_mod$data[, TISSUE := factor(TISSUE, levels=unique(data_mod$data$TISSUE))]

       # plot
       p <- ggplot(data    = data_mod$data,
                   mapping = aes(x = TISSUE, y = STUDY, fill=TPM)) +
         geom_tile() +
         viridis::scale_fill_viridis(option="inferno", direction=1) +
         labs(x = "Tissue", y = "Study", fill="TPM")


     #--------------------------------------------
     # Tissue expression (TPM of different tissues) plot - for eQTL data
     } else if(input$sensitivity_plot == "Tissue - clusters") {

       validate(
         need(nrow(unique(data_mod$data[, list(STUDY,TISSUE)]))>1, '>1 dataset required for PCA cluster analysis'),
       )

       # Calculate PCA / clustering datatable
       beta_matrix <- data_mod$data[QUANT_METHOD=="ge" & grepl("^rs[0-9]+$",RSID), list(RSID, DATASET=paste(STUDY,"-",TISSUE), BETA)]
       beta_matrix <- beta_matrix[, .SD[which.max(BETA)], by=c("RSID","DATASET")] # a few duplicates in a dataset, remove
       beta_matrix <- tidyr::complete(beta_matrix, RSID, DATASET) |> data.table::as.data.table()
       beta_matrix <- data.table::dcast(beta_matrix, RSID ~ DATASET, value.var="BETA")
       beta_matrix <- stats::na.omit(beta_matrix)
       t_beta_matrix <- data.table::as.data.table( t(beta_matrix[,-1]) )
       data.table::setnames(t_beta_matrix, names(t_beta_matrix), beta_matrix$RSID)
       t_beta_matrix[, DATASET := names(beta_matrix)[-1]]
       data.table::setcolorder(t_beta_matrix, "DATASET")

       # run PCA
       pca <- stats::prcomp(t_beta_matrix[,-1], center=TRUE, scale=TRUE)

       # plot
       p <- autoplot(pca, data=t_beta_matrix, color='DATASET') +
         geom_point(aes(color=DATASET), size=5) +
         theme_classic() +
         labs(color="Dataset")

     }


      # return
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
