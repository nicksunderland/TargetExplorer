#' basic_coloc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets
#' @importFrom shinyalert shinyalert
#'
mod_basic_coloc_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6, prettyRadioButtons(inputId=ns("point_select1"), label="Select", choices=c("all","index", "group"), selected="index", inline=FALSE)),
      column(6, prettyRadioButtons(inputId=ns("point_select2"), label="Select", choices=c("all","index", "group"), selected="index", inline=FALSE)),
    ),
    fluidRow(
      column(12, prettyRadioButtons(inputId=ns("beta_dir"), label="\u03B2 concordance", choices=c("any","concordant", "disconcordant"), selected="any", inline=TRUE)),
    )
  ) # tagList end
}


#' r_susier Server Functions
#' @description
#' The aim of this module is to return a data.table into `data_module$data` with the
#' following structure: [RSID_1, RSID_2, BP_1, BP_2, nlog10P_1, nlog10P_2, group, index]
#'
#' @noRd
mod_basic_coloc_server <- function(id, source_1_module, source_2_module, data_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #==========================================
    # observe selectors
    #==========================================
    session$userData[[ns("plot_dataset")]] <- observe({

      req(source_1_module$data, source_2_module$data, input$point_select1, input$point_select2, input$beta_dir)

      # harmonise the datasets
      harm <- genepi.utils::harmonise(gwas1 = source_1_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P)],
                                      gwas2 = source_2_module$data[, list(SNP=RSID,RSID,CHR,BP,EA,OA,EAF,BETA,SE,P,nlog10P)],
                                      gwas1_trait = "1",
                                      gwas2_trait = "2",
                                      merge = c("SNP"="SNP"))

      # keep only what we need for plotting
      harm <- harm[keep==TRUE, list(RSID_1,BP_1,BETA_1,nlog10P_1,RSID_2,BP_2,BETA_2,nlog10P_2)]

      # get the variants from source 1
      if(input$point_select1!="all" && !all(c("index","group") %in% names(source_1_module$data))) {
        selector <- "not_present"
      } else {
        selector <- input$point_select1
      }
      s1_variants <- switch(selector,
                            "all"         = { source_1_module$data[,              list(RSID)] },
                            "index"       = { source_1_module$data[index==TRUE,   list(RSID)] },
                            "group"       = { source_1_module$data[!is.na(group), list(RSID)] },
                            "not_present" = { source_1_module$data[rep(F,.N),     list(RSID)] })

      # get the variants from source 2
      if(input$point_select2!="all" && !all(c("index","group") %in% names(source_2_module$data))) {
        selector <- "not_present"
      } else {
        selector <- input$point_select2
      }
      s2_variants <- switch(selector,
                            "all"         = { source_2_module$data[,              list(RSID)] },
                            "index"       = { source_2_module$data[index==TRUE,   list(RSID)] },
                            "group"       = { source_2_module$data[!is.na(group), list(RSID)] },
                            "not_present" = { source_2_module$data[rep(F,.N),     list(RSID)] })

      # join
      all_variants <- data.table::merge.data.table(s1_variants, s2_variants, by="RSID", all=FALSE)
      all_variants[, c("index","group") := list(TRUE,factor(1:.N))]
      harm <- harm[all_variants, on=c("RSID_1"="RSID")]
      harm[, RSID := RSID_1] # need to maintain an RSID column, e.g. for MR module to filter on

      # assess beta concordance if needed
      if(input$beta_dir == "concordant") {

        harm <- harm[sign(BETA_1)==sign(BETA_2), ]

      } else if(input$beta_dir == "disconcordant") {

        harm <- harm[sign(BETA_1)!=sign(BETA_2), ]

      }

      # add to main data slot for plotting
      data_module$data <- NULL
      data_module$data <- harm

    })


    #==========================================
    # reset grouping
    #==========================================
    session$userData[[ns("reset")]] <- observeEvent(input$reset, {

      if(all(c("group","index") %in% names(data_module$data))) {
        data_module$data[, c("group","index") := NULL]
        tmp <- data.table::copy(data_module$data)
        data_module$data <- NULL
        data_module$data <- tmp
      }

    })




  })
}
