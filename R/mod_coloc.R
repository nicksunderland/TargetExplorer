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
                                        selectInput(inputId = ns("source_1"),
                                                    label   = "Dataset 1",
                                                    choices = c("")),
                                        prettyRadioButtons(inputId = ns("source_1_type"),
                                                           label   = "Type",
                                                           choices = c("quant","cc"),
                                                           selected = "quant",
                                                           inline   = TRUE)),
                                 column(6,
                                        selectInput(inputId = ns("source_2"),
                                                    label   = "Dataset 2",
                                                    choices = c("")),
                                        prettyRadioButtons(inputId = ns("source_2_type"),
                                                           label   = "Type",
                                                           choices = c("quant","cc"),
                                                           selected = "quant",
                                                           inline   = TRUE))
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
                                 column(6, selectInput(inputId  = ns("ld_reference"),
                                                       label    = "LD reference",
                                                       choices  = c("EUR_1000G"))),
                                 column(6, prettyRadioButtons(inputId  = ns("method"),
                                                              label    = "Method",
                                                              choices  = c("Single causal","SuSiE"),
                                                              selected = "Single causal"))
                               ),
                               fluidRow(
                                 column(6, actionButton(inputId = ns("run_finemap"),
                                                        label   = "Finemap")),
                                 column(6, actionButton(inputId = ns("run_coloc"),
                                                        label   = "Coloc"))
                               ),
                               hr(),
                               fluidRow(
                                 column(6,
                                        sliderTextInput(inputId  = ns("credible_set_p"),
                                                        label    = "Credible set %",
                                                        choices  = c(0.50,0.60,0.70,0.80,0.90,0.95,0.99,1.0),
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
    # Data module server for the Coloc module
    #==========================================
    data_mod   <- mod_data_server(id="data", gene_module=app$modules$gene)
    remove_mod <- mod_remove_server(id="remove", app=app, parent_id=id, parent_inputs=input)


    #==========================================
    # Run coloc button
    #==========================================
    session$userData[[ns("run_coloc")]] <- observeEvent(input$run_coloc, {

      # check data
      if(is.null(app$modules[[input$source_1]]$data) || is.null(app$modules[[input$source_2]]$data)) return(NULL)






      # HARMONISE!!!






      # source 1
      d1 <- app$modules[[input$source_1]]$data

      # test for and deal with duplicate IDs
      if(sum(duplicated(d1$RSID), na.rm=TRUE) > 0) {
        d1 <- d1[d1[, .I[which.min(P)], by=RSID]$V1]
        showNotification("Duplicate RSIDs found in source 1, taking the lowest P-value variant", type="warning")
      }

      # create the coloc object
      D1 <- list(
        snp      = d1$RSID,
        position = d1$BP,
        MAF      = d1$EAF,
        beta     = d1$BETA,
        varbeta  = d1$SE ^ 2,
        N        = d1$N,
        type     = input$source_1_type
      )

      # check the dataset
      tryCatch({
        coloc::check_dataset(d=D1, suffix="source 1")
      },
      error=function(e) {
        showNotification(paste0("Source 1 dataset check failed - ", e), type="error", duration=10)
        return(NULL)
      },
      warning=function(w) {
        showNotification(paste0("Source 1 dataset check warning - ", w), type="warning", duration=10)
      })

      # source 2
      d2 <- app$modules[[input$source_2]]$data

      # test for and deal with duplicate IDs
      if(sum(duplicated(d2$RSID), na.rm=TRUE) > 0) {
        d2 <- d2[d2[, .I[which.min(P)], by=RSID]$V1]
        showNotification("Duplicate RSIDs found in source 2, taking the lowest P-value variant", type="warning")
      }

      # create the coloc object
      D2 <- list(
        snp      = d2$RSID,
        position = d2$BP,
        MAF      = d2$EAF,
        beta     = d2$BETA,
        varbeta  = d2$SE ^ 2,
        N        = d2$N,
        type     = input$source_2_type
      )

      # check the dataset
      tryCatch({
        coloc::check_dataset(d=D2, suffix="source 2")
      },
      error=function(e) {
        showNotification(paste0("Source 2 dataset check failed - ", e), type="error", duration=10)
        return(NULL)
      },
      warning=function(w) {
        showNotification(paste0("Source 2 dataset check warning - ", w), type="warning", duration=10)
      })

      # run the colocalisation function and put the result in data2
      data_mod$data2 <- coloc::coloc.abf(dataset1=D1, dataset2=D2)

      # calculate the cumulative H4.posterior probability
      data_mod$data2$results <- data.table::as.data.table(data_mod$data2$results)
      data_mod$data2$results[order(SNP.PP.H4, decreasing=TRUE), cumsum_pp_h4 := cumsum(SNP.PP.H4)]

      # depending on the decision parameters add a coloc flag to which ever dataset is going to be the downstream
      if(input$downstream_dataset == "Dataset 1") {

        dat <- data.table::copy(app$modules[[input$source_1]]$data)

      # "Dataset 2"
      } else {

        dat <- data.table::copy(app$modules[[input$source_2]]$data)

      }

      # join cumsum posterior probability
      dat[data_mod$data2$results, cumsum_pp_h4 := i.cumsum_pp_h4, on=c("RSID"="snp")]

      # determine credible set - flag as coloc (i.e. when is the cumulative probability > 95%)
      dat[order(cumsum_pp_h4), coloc := ifelse(.I <= which(cumsum_pp_h4>=input$credible_set_p)[1], TRUE, FALSE)]

      # assign to data module `data`
      data_mod$data <- dat
    })


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
    session$userData[[ns("method")]] <- observeEvent(input$method, {
      if(input$method == "SuSiE") {
        shinyjs::enable("ld_reference")
      } else if(input$method == "Single causal") {
        shinyjs::disable("ld_reference")
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
    # Locus plot
    #==========================================
    output$locus_plot <- renderPlot({

      print("calling locus plot")

      # check data
      validate(
        need(!is.null(app$modules[[input$source_1]]$data) | !is.null(app$modules[[input$source_2]]$data), 'No data imported for either dataset')
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
      if(!is.null(app$modules[[input$source_1]]$data)) {

        y_max[[1]] <- max(app$modules[[input$source_1]]$data$nlog10P, na.rm=TRUE)

        p <- p +
          annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05)

        # add downstream label if UI indicates
        if(input$downstream_dataset=="Dataset 1") {
          p <- p +
            annotate(geom = "text", x=app$modules$gene$end+(app$modules$gene$flanks_kb*1000)-1e4, y=ceiling(y_max[[1]]), label="downstream", color="darkred", alpha=0.5)
        }

        # see if there are variables to colour the points by
        if(all(c("QUANT_METHOD","STUDY_TISSUE") %in% names(app$modules[[input$source_1]]$data))) {

          p <- p +
            geom_point(data    = app$modules[[input$source_1]]$data,
                       mapping = aes(x=BP, y=nlog10P, color=STUDY_TISSUE, shape=QUANT_METHOD)) +
            labs(color="Tissue", shape="Quant. method")

        } else {

          p <- p +
            geom_point(data    = app$modules[[input$source_1]]$data,
                       mapping = aes(x=BP, y=nlog10P),
                       color   = "lightgray")

        }

      }

      # plot source 2 data negative
      if(!is.null(app$modules[[input$source_2]]$data)) {

        y_max[[2]] <- max(app$modules[[input$source_2]]$data$nlog10P, na.rm=TRUE)

        p <- p +
          annotate(geom = "rect", xmin=app$modules$gene$start, xmax=app$modules$gene$end, ymin=-Inf, ymax=0, fill="blue", alpha = 0.05)

        # add downstream label if UI indicates
        if(input$downstream_dataset=="Dataset 2") {
          p <- p +
            annotate(geom = "text", x=app$modules$gene$end+(app$modules$gene$flanks_kb*1000)-1e4, y=-ceiling(y_max[[2]]), label="downstream", color="darkred", alpha=0.5)
        }

        # see if there are variables to colour the points by
        if(all(c("QUANT_METHOD","STUDY_TISSUE") %in% names(app$modules[[input$source_2]]$data))) {

          p <- p +
            geom_point(data    = app$modules[[input$source_2]]$data,
                       mapping = aes(x=BP, y=-nlog10P, color=STUDY_TISSUE, shape=QUANT_METHOD)) +
            labs(color="Tissue", shape="Quant. method")

        } else {

          p <- p +
            geom_point(data    = app$modules[[input$source_2]]$data,
                       mapping = aes(x=BP, y=-nlog10P),
                       color   = "darkgray")

        }

      }

      # if two sets of data draw the x-axis at y=0
      if(!is.null(app$modules[[input$source_1]]$data) && !is.null(app$modules[[input$source_2]]$data)) {
        p <- p + geom_hline(yintercept = 0)
      }

      # see if there is colocalised data to plot
      if(!is.null(data_mod$data) && "cumsum_pp_h4" %in% names(data_mod$data)) {

        # recalculate the credible set
        data_mod$data[order(cumsum_pp_h4), coloc := ifelse(.I <= which(cumsum_pp_h4>input$credible_set_p)[1], TRUE, FALSE)]

        # if source 1 data - add the coloc points
        if(!is.null(app$modules[[input$source_1]]$data)) {

          p <- p +
            geom_point(data    = app$modules[[input$source_1]]$data[RSID %in% data_mod$data[coloc==TRUE, RSID], ],
                       mapping = aes(x=BP, y=nlog10P),
                       color   = "red", fill = "red", shape=24, size=3) +
            geom_label_repel(data    = app$modules[[input$source_1]]$data[RSID %in% data_mod$data[coloc==TRUE, RSID], ],
                             mapping = aes(label=RSID,  x=BP, y=nlog10P),
                             max.overlaps = Inf)

        }

        # if source 2 data - add the coloc points
        if(!is.null(app$modules[[input$source_2]]$data)) {

          p <- p +
            geom_point(data    = app$modules[[input$source_2]]$data[RSID %in% data_mod$data[coloc==TRUE, RSID], ],
                       mapping = aes(x=BP, y=-nlog10P),
                       color   = "red", fill = "red", shape=24, size=3) +
            geom_label_repel(data    = app$modules[[input$source_2]]$data[RSID %in% data_mod$data[coloc==TRUE, RSID], ],
                             mapping = aes(label=RSID,  x=BP, y=-nlog10P),
                             max.overlaps = Inf)

        }

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
