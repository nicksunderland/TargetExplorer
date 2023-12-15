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
                  sidebarPanel(p(strong(paste0("Controls [",id,"]"))),
                               width = 3,
                               hr(),
                               fluidRow(
                                 column(6, selectInput(inputId = ns("source_1"),
                                                       label   = "Dataset 1",
                                                       choices = c(""))),
                                 column(6, selectInput(inputId = ns("source_2"),
                                                       label   = "Dataset 2",
                                                       choices = c("")))
                               ),
                               fluidRow(
                                 column(6, prettyCheckboxGroup(inputId = ns("mr_method"),
                                                               label   = "MR method",
                                                               choices = c("mr_wald_ratio","mr_egger_regression","mr_weighted_median",
                                                                            "mr_ivw","mr_simple_mode", "mr_weighted_mode"),
                                                               selected= c("mr_wald_ratio", "mr_egger_regression", "mr_weighted_median",
                                                                            "mr_ivw", "mr_simple_mode", "mr_weighted_mode"),
                                                               shape   = "curve",
                                                               inline  = FALSE)),
                                 column(6, actionButton(inputId = ns("run_mr"),
                                                        label   = "Run MR"))
                               ),
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            p(strong("Mendelian Randomisation:")),
                            plotOutput(outputId = ns("mr_plot"),
                                       height   = "500px",
                                       brush    = ns("mr_plot_brush")),
                            tableOutput(outputId = ns("mr_plot_table")),
                  ) # main panel end
    ), # sidebar layout end
    hr()
  ) # div end
}


#' MR Server Functions
#' @import ggplot2 ggrepel
#' @noRd
mod_mr_server <- function(id, app){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # R CMD checks
    # BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- NULL

    #==========================================
    # Data module server for the GWAS module
    #==========================================
    data <- mod_data_server(id="data", app$modules$gene)

    data$data <- data.frame(x=runif(1000), y=runif(1000), color=seq(0,1,length.out=100))


    #==========================================
    # Manhattan / LocusZoom plot
    #==========================================
    output$mr_plot <- renderPlot({

      # check data
      validate(
        need(!is.null(data$data), 'No MR results data currently')
      )

      # create the locus plot
      p <- ggplot(data    = data$data,
                  mapping = aes(x=x, y=y, color=color)) +
        geom_point() +
        # theme_classic() +
        # lims(x = c(gene$start - gene$flanks_kb*1000, gene$end + gene$flanks_kb*1000)) +
        labs(x        = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y        = expression(paste("-log"[10], plain(P))))

      return(p)
    })


    #==========================================
    # Point select table - for brushed point
    #==========================================
    output$mr_plot_table <- renderTable({
      tryCatch({
        bp <- brushedPoints(data$data, input$mr_plot_brush)
        if(nrow(bp)>0) return(bp) else return(NULL)
      },
      error = function(e) {
        return(NULL)
      })
    })


    #==========================================
    # Run MR button
    #==========================================
    observeEvent(input$run_mr, {

      data$data[sample(1:100, size=25), "y"] <- data$data[sample(1:100, size=25), "y"]*runif(25)
    #   # check data
    #   if(is.null(data$data)) return(NULL)
    #
    #   # get the reference file
    #   plink_ref <- make_1000G_ref_subset(chrom = gene$chr,
    #                                      from  = gene$start - gene$flanks_kb*1000,
    #                                      to    = gene$end + gene$flanks_kb*1000)
    #
    #
    #   # reset if previously run (issues with factors being reset)
    #   if("clump" %in% names(data$data)) {
    #
    #     data$data$clump <- NULL
    #     data$data$index <- NULL
    #
    #   }
    #
    #   # run clumping
    #   data$data <- genepi.utils::clump(gwas      = data$data,
    #                            p1        = input$clump_p1,
    #                            p2        = input$clump_p2,
    #                            r2        = input$clump_r2,
    #                            kb        = input$clump_kb,
    #                            plink2    = get_plink2_exe(),
    #                            plink_ref = plink_ref) |> as.data.frame()
    #
    })



  })
}
