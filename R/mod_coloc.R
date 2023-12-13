#' #' step UI Function
#' #'
#' #' @description A shiny Module.
#' #'
#' #' @param id,input,output,session Internal parameters for {shiny}.
#' #'
#' #' @noRd
#' #' @import shiny shinyWidgets
#' #'
#' mod_gwas_ui <- function(id){
#'   ns <- NS(id)
#'   # need to wrap in div and give an id in order to remove (https://www.youtube.com/watch?app=desktop&v=W7ES6QYvN_c)
#'   div(
#'     id = id,
#'     sidebarLayout(position = "right",
#'                   sidebarPanel(p(strong("Controls")),
#'                                width = 3,
#'                                hr(),
#'                                fluidRow(
#'                                  column(6, selectInput(inputId  = ns("data_source"),
#'                                                        label    = "Data source",
#'                                                        choices  = c("", "file input"),
#'                                                        selected = "")),
#'                                  column(6, fileInput(inputId = ns("gwas_file"),
#'                                                      label   = "File",
#'                                                      buttonLabel = "file path...")),
#'                                ),
#'                                fluidRow(
#'                                  column(8, prettyRadioButtons(inputId = ns("gwas_type"),
#'                                                               label   = "Type",
#'                                                               choices = c("Exposure", "Outcome"),
#'                                                               inline  = TRUE)),
#'                                  column(4, actionButton(inputId = ns("import"),
#'                                                         label   = "Import")),
#'                                  tags$style(type='text/css', paste0("#",ns("import")," { width:100%; margin-top: 25px;}"))
#'                                ),
#'                                hr(),
#'                                fluidRow(
#'                                  column(6, p(strong("Clumping"))),
#'                                  column(6, actionButton(inputId = ns("clump"),
#'                                                         label   = "Clump data")),
#'                                ),
#'                                fluidRow(
#'                                  column(6, sliderTextInput(inputId  = ns("clump_p1"),
#'                                                            label    = "p1",
#'                                                            choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
#'                                                            selected = 5e-8,
#'                                                            grid     = TRUE)),
#'                                  column(6, sliderTextInput(inputId  = ns("clump_p2"),
#'                                                            label    = "p2",
#'                                                            choices  = c(5e-8,5e-6,5e-4,0.01,0.05,0.1,0.5,1.0),
#'                                                            selected = 1.0,
#'                                                            grid     = TRUE))
#'                                ),
#'                                fluidRow(
#'                                  column(6, sliderTextInput(inputId  = ns("clump_r2"),
#'                                                            label    = "r2",
#'                                                            choices  = c(0.0001,0.001,0.01,seq(0.1,1,0.1)),
#'                                                            selected = 0.001,
#'                                                            grid     = TRUE)),
#'                                  column(6, sliderTextInput(inputId  = ns("clump_kb"),
#'                                                            label    = "kb",
#'                                                            choices  = seq(0,1000,50),
#'                                                            selected = 250,
#'                                                            grid     = TRUE))
#'                                ),
#'                   ), # sidebar panel end
#'                   mainPanel(width = 9,
#'                             # Locus zoom plot
#'                             p(strong("Locus plot:")),
#'                             plotOutput(outputId = ns("locus_plot"),
#'                                        height   = "500px",
#'                                        brush    = ns("locus_plot_brush")),
#'                             tableOutput(outputId = ns("locus_plot_table")),
#'                   ) # main panel end
#'     ), # sidebar layout end
#'     hr()
#'   ) # div end
#' }
#'
#'
#' #' gwas Server Functions
#' #' @import ggplot2 ggrepel
#' #' @noRd
#' mod_gwas_server <- function(id, base_module){
#'   moduleServer( id, function(input, output, session){
#'     ns <- session$ns
#'
#'
#'     # reactive values for the GWAS module
#'     v <- reactiveValues(
#'       genes = NULL
#'     )
#'
#'
#'     # MAIN LOCUS OUTPUT PLOT
#'     output$locus_plot <- renderPlot({
#'
#'       # check data
#'       validate(
#'         need(!is.null(base_module$data), 'No data imported, click the import button')
#'       )
#'
#'       # colours
#'       color_list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#'
#'       # create the locus plot
#'       p <- ggplot(data    = base_module$data,
#'                   mapping = aes(x=BP, y=log10P)) +
#'         geom_point(color="lightgray") +
#'         annotate(geom = "rect", xmin=base_module$gene_start, xmax=base_module$gene_end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
#'         theme_classic() +
#'         lims(x = c(base_module$gene_start - base_module$gene_flanks_kb*1000, base_module$gene_end + base_module$gene_flanks_kb*1000)) +
#'         labs(x        = paste0("Chromosome ", base_module$gene_chr, " position"),
#'              y        = expression(paste("-log"[10], plain(P))),
#'              subtitle = paste0(input$gwas_type, ": ", input$data_source))
#'
#'       # if gene data then add
#'       if(!is.null(v$genes)) {
#'         p <- p +
#'           geom_hline(yintercept = 0) +
#'           geom_rect(data = v$genes[v$genes$STRAND=="+",  ],
#'                     inherit.aes=FALSE,
#'                     mapping = aes(xmin=BP_START, xmax=BP_END, ymin=-1.75, ymax=-1.25), fill="grey", alpha=0.3) +
#'           geom_rect(data = v$genes[v$genes$STRAND=="-",  ],
#'                     inherit.aes=FALSE,
#'                     mapping = aes(xmin=BP_START, xmax=BP_END, ymin=-2.5, ymax=-2), fill="grey", alpha=0.3) +
#'           geom_text_repel(data = v$genes[v$genes$STRAND=="-",  ],
#'                           mapping = aes(label = GENE_NAME, x=(BP_END-BP_START)/2 + BP_START, y =-2.25),
#'                           direction = "x", min.segment.length = 0.25) +
#'           geom_text_repel(data = v$genes[v$genes$STRAND=="+",  ],
#'                           mapping = aes(label = GENE_NAME, x=(BP_END-BP_START)/2 + BP_START, y =-1.5),
#'                           direction = "x", min.segment.length = 0.25)
#'       }
#'
#'       # if there is clumped data, plot
#'       if(all(c("index","clump") %in% colnames(base_module$data))) {
#'         p <- p +
#'           geom_point(data = base_module$data[!is.na(base_module$data$clump), ],            mapping = aes(x=BP, y=log10P, color=clump, fill=clump), shape=23) +
#'           geom_vline(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(xintercept=BP), linetype="dotted", color="darkred") +
#'           geom_point(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(x=BP, y=log10P), size=3, fill="red", color="red", shape=24) +
#'           geom_label(data = base_module$data[which(base_module$data$index==TRUE), ],       mapping = aes(label=clump, x=BP, y=-0.5)) +
#'           geom_label_repel(data = base_module$data[which(base_module$data$index==TRUE), ], mapping = aes(label=RSID,  x=BP, y=log10P)) +
#'           labs(color = "Clump", fill = "Clump")
#'       }
#'
#'       return(p)
#'     })
#'
#'
#'     # Table for brushed (selected) locus plot points
#'     output$locus_plot_table <- renderTable({
#'       tryCatch({
#'         bp <- brushedPoints(base_module$data, input$locus_plot_brush)
#'         if(nrow(bp)==0 && !is.null(v$genes)) {
#'           bp <- v$genes[(input$locus_plot_brush$xmin < v$genes$BP_START) & (input$locus_plot_brush$xmax > v$genes$BP_END) &
#'                           ((input$locus_plot_brush$ymin < -1.75 & input$locus_plot_brush$ymax > -1.25 & v$genes$STRAND=="+") |
#'                              (input$locus_plot_brush$ymin < -2.50 & input$locus_plot_brush$ymax > -2.00 & v$genes$STRAND=="-")), ]
#'         }
#'         if(nrow(bp)>0) return(bp) else return(NULL)
#'       },
#'       error = function(e) {
#'         return(NULL)
#'       })
#'     })
#'
#'
#'     # observe the import button
#'     observeEvent(input$import, {
#'
#'       # nothing selected
#'       if(input$data_source=="") return()
#'
#'       # custom file input
#'       if(input$data_source=="file input") {
#'
#'         if(is.null(input$gwas_file)) {
#'           return()
#'         } else {
#'           base_module$data <- input$gwas_file
#'           base_module$data$log10P <- -log10(base_module$data$P)
#'         }
#'
#'         # internal file input
#'       } else {
#'
#'         # GWAS data lives in the base module
#'         base_module$data <- read_internal_data(type   = "gwas",
#'                                                source = input$data_source,
#'                                                chrom  = base_module$gene_chr,
#'                                                start  = base_module$gene_start - base_module$gene_flanks_kb*1000,
#'                                                end    = base_module$gene_end + base_module$gene_flanks_kb*1000)
#'         base_module$data$log10P <- -log10(base_module$data$P)
#'
#'         # Gene info lives just in the GWAS module
#'         v$genes <- read_gene_data(source = "gencode",
#'                                   chrom  = base_module$gene_chr,
#'                                   start  = base_module$gene_start - base_module$gene_flanks_kb*1000,
#'                                   end    = base_module$gene_end + base_module$gene_flanks_kb*1000)
#'
#'       }
#'
#'     })
#'
#'
#'     # observe the clump button
#'     observeEvent(input$clump, {
#'
#'       # check data
#'       if(is.null(base_module$data)) return(NULL)
#'
#'       # get the reference file
#'       plink_ref <- make_1000G_ref_subset(chrom = base_module$gene_chr,
#'                                          from  = base_module$gene_start - base_module$gene_flanks_kb*1000,
#'                                          to    = base_module$gene_end + base_module$gene_flanks_kb*1000)
#'
#'
#'       # reset if previously run (issues with factors being reset)
#'       if("clump" %in% names(base_module$data)) {
#'
#'         base_module$data$clump <- NULL
#'         base_module$data$index <- NULL
#'
#'       }
#'
#'       # run clumping
#'       base_module$data <- genepi.utils::clump(gwas      = base_module$data,
#'                                               p1        = input$clump_p1,
#'                                               p2        = input$clump_p2,
#'                                               r2        = input$clump_r2,
#'                                               kb        = input$clump_kb,
#'                                               plink2    = get_plink2_exe(),
#'                                               plink_ref = plink_ref) |> as.data.frame()
#'
#'     })
#'
#'
#'     # observe the data source select input
#'     observeEvent(input$data_source, {
#'
#'       if(input$data_source=="") {
#'         # update at the beginning (when the choice is "")
#'         files <- get_available_data_sources("gwas")
#'         updateSelectInput(session, inputId="data_source", choices=c("hba1c_jurgens2022", names(files)))
#'       } else if (input$data_source=="file input") {
#'         # using user dataset, enable file input
#'         shinyjs::enable(id="gwas_file")
#'       } else {
#'         # using example dataset, disable file input
#'         shinyjs::disable(id="gwas_file")
#'       }
#'
#'     })
#'
#'
#'
#'   })
#' }
