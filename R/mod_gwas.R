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
                  sidebarPanel(p(strong(paste0("Controls [",id,"]"))),
                               width = 3,
                               hr(),
                               mod_data_ui(id=paste0(id,"-data")),
                               hr(),
                               mod_clump_ui(id=ns("clump"))
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            p(strong("Locus plot:")),
                            plotOutput(outputId = ns("locus_plot"),
                                       height   = "500px",
                                       brush    = ns("locus_plot_brush")),
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
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL

    #==========================================
    # Data and clump module servers for the GWAS module
    #==========================================
    data_mod  <- mod_data_server(id="data", gene_module=app$modules$gene)
    clump_mod <- mod_clump_server(id="clump", gene_module=app$modules$gene, data_module=data_mod)

    #==========================================
    # Manhattan / LocusZoom plot
    #==========================================
    output$locus_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(data_mod$data), 'No data imported, click the import button')
      )

      # colours
      color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

      # create the locus plot
      p <- ggplot(data    = data_mod$data,
                  mapping = aes(x=BP, y=nlog10P)) +
        geom_hline(yintercept=7.3, linetype="dotted", color="lightgrey") +
        geom_point(color="lightgray") +
        annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        lims(x = c(app$modules$gene$start - app$modules$gene$flanks_kb*1000, app$modules$gene$end + app$modules$gene$flanks_kb*1000),
             y = c(-2.5, max(8.0, max(data_mod$data$nlog10P)))) +
        labs(x        = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y        = expression(paste("-log"[10], plain(P))))

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

      # if there is clumped data, plot
      if(all(c("index","clump") %in% colnames(data_mod$data))) {
        p <- p +
          geom_point(data = data_mod$data[!is.na(data_mod$data$clump), ],            mapping = aes(x=BP, y=nlog10P, color=clump, fill=clump), shape=23) +
          geom_vline(data = data_mod$data[which(data_mod$data$index==TRUE), ],       mapping = aes(xintercept=BP), linetype="dotted", color="darkred") +
          geom_point(data = data_mod$data[which(data_mod$data$index==TRUE), ],       mapping = aes(x=BP, y=nlog10P), size=3, fill="red", color="red", shape=24) +
          geom_label(data = data_mod$data[which(data_mod$data$index==TRUE), ],       mapping = aes(label=clump, x=BP, y=-0.5)) +
          geom_label_repel(data = data_mod$data[which(data_mod$data$index==TRUE), ], mapping = aes(label=RSID,  x=BP, y=nlog10P)) +
          labs(color = "Clump", fill = "Clump")

      # scenario when 'all' checkbox is clicked - no clump column, but a TRUE index column
      } else if("index" %in% colnames(data_mod$data)) {
        p <- p +
          geom_point(data = data_mod$data,
                     mapping = aes(x=BP, y=nlog10P,
                                   color=factor(index, levels=c(TRUE), labels=c("Use all")),
                                   fill =factor(index, levels=c(TRUE), labels=c("Use all"))), size=1, shape=24) +
          scale_fill_manual(values=c("red")) +
          scale_color_manual(values=c("red")) +
          guides(color = guide_legend(title = NULL),
                 fill  = guide_legend(title = NULL))
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
    # reset clumping
    #==========================================
    observeEvent(input$reset, {

      if(all(c("clump","index") %in% names(data_mod$data))) {
        data_mod$data[, c("clump","index") := NULL]
        tmp <- data.table::copy(data_mod$data)
        data_mod$data <- NULL
        data_mod$data <- tmp
      }

    })


    #==========================================
    # Clump button
    #==========================================
    observeEvent(input$clump, {

      # check data
      if(is.null(data_mod$data)) return(NULL)

      # clean up previous runs if data present
      if("clump" %in% names(data_mod$data)) {
        data_mod$data$clump <- NULL
      }
      if("index" %in% names(data_mod$data)) {
        data_mod$data$index <- NULL
      }

      # get the reference file
      plink_ref <- make_1000G_ref_subset(chrom = app$modules$gene$chr,
                                         from  = app$modules$gene$start - app$modules$gene$flanks_kb*1000,
                                         to    = app$modules$gene$end + app$modules$gene$flanks_kb*1000)


      # reset if previously run (issues with factors being reset)
      if("clump" %in% names(data_mod$data)) {

        data_mod$data$clump <- NULL
        data_mod$data$index <- NULL

      }

      shiny::withProgress(message = 'Clumping data', value = 0, {

        shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))

        # run clumping
        data_mod$data <- genepi.utils::clump(gwas      = data_mod$data,
                                             p1        = input$clump_p1,
                                             p2        = input$clump_p2,
                                             r2        = input$clump_r2,
                                             kb        = input$clump_kb,
                                             plink2    = get_plink2_exe(),
                                             plink_ref = plink_ref) #|> as.data.frame()

        if(any(data_mod$data$index)) {
          shiny::incProgress(3/4, detail = paste("Complete"))
        } else {
          shiny::incProgress(2/4, detail = paste("Failed"))
          Sys.sleep(1)
          shiny::incProgress(1/4, detail = paste("Failed"))
        }

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
