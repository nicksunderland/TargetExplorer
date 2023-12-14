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
                  sidebarPanel(p(strong("Controls")),
                               width = 3,
                               hr(),
                               fluidRow(
                                 column(6, selectInput(inputId  = ns("data_source"),
                                                       label    = "Data source",
                                                       choices  = c("", "file input", "Local", "EBI eQTL Catalogue"),
                                                       selected = ""),
                                           fileInput(inputId = ns("gwas_file"),
                                                     label   = "File",
                                                     buttonLabel = "file path..."),
                                           prettyRadioButtons(inputId = ns("build"),
                                                              label   = "Build",
                                                              choices = c("GRCh37","GRCh38"),
                                                              selected= "GRCh37",
                                                              inline  = TRUE),
                                           actionButton(inputId = ns("import"),
                                                        label   = "Import"),
                                        ),
                                 column(6,
                                        datamods::select_group_ui(id = ns("datasets_filter"),
                                                                  params = list(
                                                                    study_label = list(inputId = "study_label", label = "Study:"),
                                                                    tissue_label= list(inputId = "tissue_label", label = "Tissue:")
                                                                  ),
                                                                  inline = FALSE),
                                        sliderTextInput(inputId  = ns("eqtl_pval"),
                                                       label    = "p-value thresh",
                                                       choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
                                                       selected = 5e-4,
                                                       grid     = TRUE)
                                        ),
                               ),
                               hr(),
                               fluidRow(
                                 column(6, p(strong("Parameters"))),
                                 column(6, actionButton(inputId = ns("run"),
                                                        label   = "Run")),
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
mod_coloc_server <- function(id, base_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # R CMD checks
    BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- nlog10p <- position <- tissue_label <- NULL

    # reactive values for the Coloc module
    v <- reactiveValues(
      ebi_eqtl_info = NULL,
      qtl_data = NULL,
      test1 = NULL,
      coloc_result = NULL,
    )

    # create the EBI study information filtering module
    ebi_dataset_module <- datamods::select_group_server(id="datasets_filter",
                                                        data = reactive(v$ebi_eqtl_info),
                                                        vars = reactive(c("study_label","tissue_label")))


    # observe the data source select input
    observeEvent(input$data_source, {

      if(input$data_source=="") {
        # update at the beginning (when the choice is "")

      } else if (input$data_source=="EBI eQTL Catalogue") {

        # query the available datasets from the API
        if(is.null(v$ebi_eqtl_info)) {

          v$ebi_eqtl_info <- get_ebi_study_info()

        }

      } else if (input$data_source=="file input") {
        # using user dataset, enable file input
        shinyjs::enable(id="gwas_file")
      } else {
        # using example dataset, disable file input
        shinyjs::disable(id="gwas_file")
      }

    })



    # observe the import button
    observeEvent(input$import, {

      if(input$data_source=="EBI eQTL Catalogue") {

        # query the API for the data
        v$qtl_data <- import_ebi_eqtl(study_info_table = v$ebi_eqtl_info,
                                      studies  = input$`datasets_filter-study_label`,
                                      tissues  = input$`datasets_filter-tissue_label`,
                                      chr      = base_module$gene_chr,
                                      bp_start = base_module$gene_start - base_module$gene_flanks_kb*1000,
                                      bp_end   = base_module$gene_end + base_module$gene_flanks_kb*1000,
                                      gene_id  = base_module$gene_id,
                                      nlog10p  = -log10(input$eqtl_pval),
                                      build    = input$build)

      }





    })


    # observe the run button
    observeEvent(input$run, {

      # library(coloc)
      # data(coloc_test_data)
      # attach(coloc_test_data)
      # v$qtl_data <- D1
      # v$test1 <- D2

      v$coloc_result <- coloc::coloc.abf(dataset1=v$qtl_data, dataset2=v$test1)

    })


    # COLOC OUTPUT TABLE
    output$coloc_plot_table <- renderTable({

      if(is.null(v$coloc_result)) return(NULL)

      t(as.data.frame(v$coloc_result$summary))

    })


    # TRAIT 1 LOCUS OUTPUT PLOT
    output$locus_plot1 <- renderPlot({

      # check data
      validate(
        need(!is.null(base_module$data), 'No data imported, click the import button')
      )

      # colours
      color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

      # create the locus plot
      p <- ggplot(data    = base_module$data,
                  mapping = aes(x=BP, y=log10P)) +
        geom_point(color="lightgray") +
        annotate(geom = "rect", xmin=base_module$gene_start, xmax=base_module$gene_end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        lims(x = c(base_module$gene_start - base_module$gene_flanks_kb*1000, base_module$gene_end + base_module$gene_flanks_kb*1000)) +
        labs(x        = paste0("Chromosome ", base_module$gene_chr, " position"),
             y        = expression(paste("-log"[10], plain(P))),
             subtitle = paste0(input$gwas_type, ": ", input$data_source))

      # if there is clumped data, plot
      if("coloc" %in% colnames(base_module$data)) {
        # p <- p +
        #   geom_point(data = base_module$data[!is.na(base_module$data$clump), ],            mapping = aes(x=BP, y=log10P, color=clump, fill=clump), shape=23) +
        #   geom_vline(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(xintercept=BP), linetype="dotted", color="darkred") +
        #   geom_point(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(x=BP, y=log10P), size=3, fill="red", color="red", shape=24) +
        #   geom_label(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(label=clump, x=BP, y=-0.5)) +
        #   geom_label_repel(data = base_module$data[which(base_module$data$index==TRUE), ], mapping = aes(label=RSID,  x=BP, y=log10P)) +
        #   labs(color = "Clump", fill = "Clump")
      }

      return(p)
    })


    # TRAIT 2 LOCUS OUTPUT PLOT
    output$locus_plot2 <- renderPlot({

      # check data
      validate(
        need(!is.null(v$qtl_data), 'No data imported, click the import button')
      )

      # colours
      color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

      # create the locus plot
      p <- ggplot(data    = v$qtl_data,
                  mapping = aes(x=position, y=nlog10p, color=tissue_label)) +
        geom_point() + #color="lightgray"
        annotate(geom = "rect", xmin=base_module$gene_start, xmax=base_module$gene_end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
        theme_classic() +
        lims(x = c(base_module$gene_start - base_module$gene_flanks_kb*1000, base_module$gene_end + base_module$gene_flanks_kb*1000)) +
        labs(x        = paste0("Chromosome ", base_module$gene_chr, " position"),
             y        = expression(paste("-log"[10], plain(P))),
             subtitle = "eQTL data")

      # if there is clumped data, plot
      if("coloc" %in% colnames(base_module$data)) {
        # p <- p +
        #   geom_point(data = base_module$data[!is.na(base_module$data$clump), ],            mapping = aes(x=BP, y=log10P, color=clump, fill=clump), shape=23) +
        #   geom_vline(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(xintercept=BP), linetype="dotted", color="darkred") +
        #   geom_point(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(x=BP, y=log10P), size=3, fill="red", color="red", shape=24) +
        #   geom_label(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(label=clump, x=BP, y=-0.5)) +
        #   geom_label_repel(data = base_module$data[which(base_module$data$index==TRUE), ], mapping = aes(label=RSID,  x=BP, y=log10P)) +
        #   labs(color = "Clump", fill = "Clump")
      }

      return(p)
    })


    # PRIOR PROBABILITY OUTPUT PLOT
    output$prob_plot1 <- renderPlot({

      # check data
      validate(
        need(!is.null(v$coloc_result), 'No data imported, click the run button')
      )

      p <- genepi.utils::plot_coloc_probabilities(v$coloc_result, rule="H4 > 0.5", type="prior")

      return(p)
    })


    # PRIOR PROBABILITY OUTPUT PLOT
    output$prob_plot2 <- renderPlot({

      # check data
      validate(
        need(!is.null(v$coloc_result), 'No data imported, click the import button')
      )

      p <- genepi.utils::plot_coloc_probabilities(v$coloc_result, rule="H4 > 0.5", type="posterior")

      return(p)
    })





  })
}
