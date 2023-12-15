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
                  sidebarPanel(p(strong(paste0("Controls [",id,"]"))),
                               width = 3,
                               hr(),
                               mod_data_ui(id=ns("data")),
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
                                 column(6,
                                        sliderTextInput(inputId  = ns("coloc_p1"),
                                                        label    = "p1",
                                                        choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                                                        selected = 5e-4,
                                                        grid     = TRUE),
                                        sliderTextInput(inputId  = ns("coloc_p12"),
                                                        label    = "p12",
                                                        choices  = c(0.0001,0.001,0.01,seq(0.1,1,0.1)),
                                                        selected = 0.001,
                                                        grid     = TRUE)),
                                 column(6,
                                        sliderTextInput(inputId  = ns("coloc_p2"),
                                                        label    = "p2",
                                                        choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                                                        selected = 5e-4,
                                                        grid     = TRUE),
                                        sliderTextInput(inputId  = ns("coloc_h4"),
                                                        label    = "H4 decision",
                                                        choices  = seq(0,1,by=0.1),
                                                        selected = 0.5,
                                                        grid     = TRUE),
                                        )
                               ),
                               fluidRow(
                                 column(6, actionButton(inputId = ns("run_finemap"),
                                                        label   = "Finemap")),
                                 column(6, actionButton(inputId = ns("run_coloc"),
                                                        label   = "Coloc")),
                               )
                  ), # sidebar panel end
                  mainPanel(width = 9,
                            # Locus zoom plot
                            p(strong("Colocalisation:")),
                            fluidRow(
                              column(6,
                                     plotOutput(outputId = ns("locus_plot1"), height = "250px"),
                                     plotOutput(outputId = ns("locus_plot2"), height = "250px")
                              ),
                              column(6,
                                     plotOutput(outputId = ns("prob_plot1"), height = "250px"),
                                     plotOutput(outputId = ns("prob_plot2"), height = "250px")
                              )
                            ),
                            tableOutput(outputId = ns("coloc_plot_table")),
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

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- nlog10P <- position <- tissue_label <- NULL

    #==========================================
    # Data module server for the Coloc module
    #==========================================
    data_mod <- mod_data_server(id="data", app$modules$gene)

    #==========================================
    # Run coloc button
    #==========================================
    observeEvent(input$run_coloc, {

      library(coloc)
      data(coloc_test_data)
      attach(coloc_test_data)

      data_mod$data2 <- coloc::coloc.abf(dataset1=D1, dataset2=D2)

    })


    #==========================================
    # Observe additions / deletions of modules
    #==========================================
    observeEvent(app$modules, {
      print("app$modules")
      print(names(app$modules))
      updateSelectInput(session, inputId="source_1", choices=names(app$modules)[names(app$modules)!="gene"])
      updateSelectInput(session, inputId="source_2", choices=names(app$modules)[names(app$modules)!="gene"])
    })

    observeEvent(list(input$source_1,input$source_2), {
      if(is.null(input$source_1) || input$source_1=="") {
        updateSelectInput(session, inputId="source_1", choices=names(app$modules)[names(app$modules)!="gene"])
      }
      if(is.null(input$source_2) || input$source_2=="") {
        updateSelectInput(session, inputId="source_2", choices=names(app$modules)[names(app$modules)!="gene"])
      }
    })


    #==========================================
    # Results table
    #==========================================
    output$coloc_plot_table <- renderTable({

      if(is.null(app$modules[[id]]$data2)) return(NULL)
      t(as.data.frame(app$modules[[id]]$data2$summary))

    })


    #==========================================
    # Trait 1 locus plot
    #==========================================
    output$locus_plot1 <- renderPlot({

      # check data
      validate(
        need(!is.null(app$modules[[input$source_1]]$data),
             paste0('No data imported for [', input$source_1, ']'))
      )

      # colours
      color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

      # create the locus plot
      p <- ggplot(data    = app$modules[[input$source_1]]$data,
                  mapping = aes(x=BP, y=nlog10P)) +
        geom_point(color="lightgray") +
        annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        lims(x = c(app$modules$gene$start - app$modules$gene$flanks_kb*1000, app$modules$gene$end + app$modules$gene$flanks_kb*1000)) +
        labs(x        = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y        = expression(paste("-log"[10], plain(P))),
             subtitle = paste0(input$gwas_type, ": ", input$data_source))

      return(p)
    })


    #==========================================
    # Trait 2 locus plot
    #==========================================
    output$locus_plot2 <- renderPlot({

      # check data
      validate(
        need(!is.null(app$modules[[input$source_2]]$data),
             paste0('No data imported for [', input$source_2, ']'))
      )

      # colours
      color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

      # create the locus plot
      p <- ggplot(data    = app$modules[[input$source_2]]$data,
                  mapping = aes(x=BP, y=nlog10P))

      if("TISSUE" %in% names(app$modules[[input$source_2]]$data)) {
        p <- p + geom_point(aes(color=TISSUE))
      } else {
        p <- p + geom_point(color="lightgray")
      }

      p <- p +
        annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        lims(x = c(app$modules$gene$start - app$modules$gene$flanks_kb*1000, app$modules$gene$end + app$modules$gene$flanks_kb*1000)) +
        labs(x = paste0("Chromosome ", app$modules$gene$chr, " position"),
             y = expression(paste("-log"[10], plain(P))),
             subtitle = "eQTL data")

      return(p)
    })


    #==========================================
    # Prior probabilities plot
    #==========================================
    output$prob_plot1 <- renderPlot({

      # check data
      validate(
        need(!is.null(app$modules[[id]]$data2), 'No data imported, click the run button')
      )

      p <- genepi.utils::plot_coloc_probabilities(app$modules[[id]]$data2, rule="H4 > 0.5", type="prior")

      return(p)
    })


    #==========================================
    # Posterior probabilities plot
    #==========================================
    output$prob_plot2 <- renderPlot({

      # check data
      validate(
        need(!is.null(app$modules[[id]]$data2), 'No data imported, click the import button')
      )

      p <- genepi.utils::plot_coloc_probabilities(app$modules[[id]]$data2, rule="H4 > 0.5", type="posterior")

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
