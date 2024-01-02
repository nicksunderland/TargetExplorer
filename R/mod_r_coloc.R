#' r_coloc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @import shiny shinyWidgets shinyalert
#'
mod_r_coloc_ui <- function(id){
  ns <- NS(id)
  cli::cli_alert_info(paste0("initialising mod_r_coloc_ui(ns=",ns("foo"),")"))
  useShinyalert()
  tagList(
    fluidRow(
      column(6,
             selectInput(inputId  = ns("func"),
                         label    = "Function",
                         choices  = c("finemap.abf","finemap.signals"),
                         selected = "finemap.signals")),
      column(4,
             actionButton(inputId = ns("run"),
                          label   = "Run")),
      column(1,
             actionButton(inputId = ns("reset"),
                          width   = "40px",
                          label   = "",
                          icon    = icon("rotate-left"))),
      tags$style(type='text/css', paste0("#",ns("run")," { width:100%; margin-top: 25px;}")),
      tags$style(type='text/css', paste0("#",ns("reset")," { width:100%; margin-top: 25px;}")),
    ),
    fluidRow(
      column(6,
             mod_reference_ui(id=ns("reference"))),
      column(3,
             prettyRadioButtons(inputId  = ns("source_type"),
                                label    = "Type",
                                choices  = c("quant","cc"),
                                selected = "quant",
                                inline   = FALSE)),
      column(3,
             prettyRadioButtons(inputId  = ns("method"),
                                label    = "Method",
                                choices  = c("mask","cond","single"),
                                selected = "mask",
                                inline   = FALSE))
    ),
    fluidRow(
      column(3,
             numericInput(inputId  = ns("l_param"),
                          label    = "L",
                          value    = 5,
                          step     = 1)),
      column(3,
             numericInput(inputId  = ns("s_param"),
                          label    = "S",
                          value    = NA_real_,
                          step     = 0.1)),
      column(3,
             numericInput(inputId  = ns("n_param"),
                          label    = "N",
                          value    = NA_real_,
                          step     = 1)),
      column(3,
             numericInput(inputId = ns("sd_y"),
                          label = "sdY",
                          value = NA_real_,
                          step  = 0.1))
    ),
    fluidRow(
      column(6,
             sliderTextInput(inputId  = ns("finemap_p1"),
                             label    = "p1",
                             choices  = c(5e-8,5e-6,1e-4,0.01,0.05,0.1,0.5,1.0),
                             selected = 1e-4,
                             grid     = TRUE)),

      column(6,
             sliderTextInput(inputId  = ns("coverage"),
                             label    = "Coverage %",
                             choices  = c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,1.0),
                             selected = 0.95,
                             grid     = TRUE))
    )
  ) # tagList end
}


#' r_coloc Server Functions
#' @noRd
mod_r_coloc_server <- function(id, gene_module, data_module){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    cli::cli_alert_info(paste0("initialising mod_r_coloc_server(ns=",ns("foo"),")"))

    # R CMD checks
    # BP <- BP_END <- BP_START <- GENE_NAME <- RSID <- clump <- log10P <- SNP <- index <- nlog10P <- NULL
    run_flag <- reactiveVal(FALSE)

    #==========================================
    # Data, clump, and remove (sub)module servers for the GWAS module
    #==========================================
    reference_mod <- mod_reference_server(id="reference", label="Reference", gene_module=gene_module, enabled=T)


    #==========================================
    # Run button
    #==========================================
    session$userData[[ns("run")]] <- observeEvent(input$run, {

      # check data
      if(is.null(data_module$data)) return(NULL)

      # if large amount of data, confirm run
      if(nrow(data_module$data) > 2000) {

        shinyalert(
          inputId = "shinyalert",
          title = "Processing time warning",
          text = paste0(nrow(data_module$data), " data points currently loaded, this may take some time. \n",
                        "Consider narrowing the genomic window and re-loading the data.\n",
                        "",
                        "Click OK to continue processing."),
          type = "warning",
          showCancelButton = TRUE,
          confirmButtonCol = "#AEDEF4",
          animation = FALSE,
          callbackR = function(ret) { if(ret) { run_flag(TRUE) } }
        )

      } else {

        run_flag(TRUE)

      }

    })


    #==========================================
    # Run function
    #==========================================
    session$userData[[ns("main")]] <- observeEvent(req(run_flag() == TRUE), {

      # reset run flag
      isolate(run_flag(FALSE))

      # reset if previously run (issues with factors being reset)
      if(any(c("group","index") %in% names(data_module$data))) {

        data_module$data$group <- NULL
        data_module$data$index <- NULL

      }

      # functions might fail
      tryCatch({

        shiny::withProgress(message = 'Finemapping data', value = 0, {

          shiny::incProgress(1/4, detail = paste("Processing", input$`filter-dataset`, "..."))

          # finemap.abf
          if(input$func == "finemap.abf") {

            # create coloc dataset
            D1 <- make_coloc_dataset(dat  = data_module$data,
                                     type = input$source_type,
                                     sdY  = input$sd_y,
                                     N    = input$n_param)

            # if assuming only one significant variant run coloc::finemap.abf()
            results <- coloc::finemap.abf(D1,
                                          p1 = input$finemap_p1)

            # update the data with grouping
            data_module$data <- calc_credible_set(result      = results,
                                                  dat_join_to = data_module$data,
                                                  coverage    = input$coverage,
                                                  result_type = "finemap.abf")

            # store the LD structure of the region
            data_module$finemap.abf_result <- data.table(RSID = results$snp,
                                                         BP   = results$position,
                                                         PP   = results$SNP.PP)

          } else if(input$func == "finemap.signals") {

            # get the reference file
            if(grepl("1kGv3", reference_mod$ref_path)) {

              plink_ref <- make_ref_subset(ref_path = reference_mod$ref_path,
                                           chrom    = gene_module$chr,
                                           from     = gene_module$start - gene_module$flanks_kb*1000,
                                           to       = gene_module$end + gene_module$flanks_kb*1000)

              dat_ld_obj <- genepi.utils::ld_matrix(dat       = data_module$data,
                                                    method    = "r",
                                                    plink2    = get_plink2_exe(),
                                                    plink_ref = plink_ref,
                                                    ukbb_ref  = NULL)

            } else if(grepl("UKB_LD", reference_mod$ref_path)) {

              # get the UKBB LD file path (downloads if not in cache)
              ukbb_ref_dt <- genepi.utils::download_ukbb_ld(chr           = gene_module$chr,
                                                            bp_start      = gene_module$start - gene_module$flanks_kb*1000,
                                                            bp_end        = gene_module$end + gene_module$flanks_kb*1000,
                                                            ukbb_ld_cache = file.path(reference_mod$ref_path, "cache"))

              # create the LD matrix for the region variants
              dat_ld_obj <- genepi.utils::ld_matrix(dat       = data_module$data,
                                                    method    = "r",
                                                    plink2    = get_plink2_exe(),
                                                    plink_ref = NULL,
                                                    ukbb_ref  = ukbb_ref_dt$root_file)
            }

            # store the r corr matrix
            data_module$ld_matrix_r  <- dat_ld_obj[["ld_mat"]]

            # reset the r2 matrix
            data_module$ld_matrix_r2 <- NULL

            # make and check the dataset
            D1 <- make_coloc_dataset(dat  = dat_ld_obj[["dat"]],
                                     type = input$source_type,
                                     sdY  = input$sd_y,
                                     S    = input$s_param,
                                     N    = input$n_param,
                                     ld   = dat_ld_obj[["ld_mat"]])

            # run finemap.signals
            results <- coloc::finemap.signals(D         = D1,
                                              LD        = D1$LD,
                                              method    = input$method,
                                              r2thr     = 0.01,
                                              sigsnps   = NULL,
                                              pthr      = input$finemap_p1,
                                              maxhits   = ifelse(input$method=="single",1,input$l_param),
                                              return.pp = FALSE)

            # update the data with grouping
            data_module$data <- calc_credible_set(result      = results,
                                                  dat_join_to = data_module$data,
                                                  coverage    = input$coverage,
                                                  result_type = "finemap.signals")

            # store data for the kriging_rss plot
            data_module$kriging_rss         <- susieR::kriging_rss(D1$beta/(D1$varbeta ^ 0.5), D1$LD, n=D1$N)$conditional_dist
            data_module$kriging_rss$RSID    <- D1$snp
            data_module$kriging_rss$outlier <- ifelse(data_module$kriging_rss$logLR  > 2 &
                                                      abs(data_module$kriging_rss$z) > 2, TRUE, FALSE)


          }

          # warn if no groups found
          if(all(is.na(data_module$data$group))) {
            shiny::incProgress(3/4, detail = paste("Failed"))
            showNotification("Finemapping failed: \nno groups were found with the provided parameters, or coloc failed", type="error")
          } else {
            shiny::incProgress(3/4, detail = paste("Complete"))
          }

        }) # end withProgress

      },
      error=function(e) {

        showNotification(paste0(input$func, " failed - ", e), type="error", duration = 10)
        return(NULL)

      }) # end tryCatch

    }) # end Run button


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


    #==========================================
    # function / method choice and input control activation
    #==========================================
    observeEvent(list(input$func, input$source_type, data_module$data, input$method), {

      req(input$func)
      req(input$source_type)
      req(input$method)

      # common params
      if(input$source_type == "quant") {

        shinyjs::disable("s_param")
        shinyjs::enable("sd_y")

      } else if(input$source_type == "cc") {

        shinyjs::enable("s_param")
        shinyjs::disable("sd_y")

      }

      # sample size
      if(!is.null(data_module$data) && "N" %in% names(data_module$data) && is.na(input$n_param)) {

        sample_size <- max(data_module$data$N, na.rm=TRUE)
        updateNumericInput(inputId="n_param", value=sample_size)

      } else {

        sample_size <- input$n_param

      }

      # sdY
      if(!is.null(data_module$data) && all(c("SE","EAF") %in% names(data_module$data)) && !is.na(sample_size)) {

        if(input$source_type == "quant") {
          sdY_estimate <- sdY.est(data_module$data$SE ^ 2,
                                  data_module$data$EAF,
                                  sample_size)
          updateNumericInput(inputId="sd_y", value=sdY_estimate)
        } else {
          updateNumericInput(inputId="sd_y", value=NA_real_)
        }

      }

      # function requirements
      if(input$func == "finemap.abf") {

        # required controls
        shinyjs::enable("source_type")
        shinyjs::enable("finemap_p1")
        shinyjs::enable("coverage")
        shinyjs::disable("reference-reference")
        shinyjs::disable("method")
        shinyjs::disable("l_param")

      } else if(input$func == "finemap.signals") {

        # required controls
        shinyjs::enable("source_type")
        shinyjs::enable("finemap_p1")
        shinyjs::enable("coverage")
        shinyjs::enable("reference-reference")
        shinyjs::enable("method")
        shinyjs::enable("l_param")

        if(input$method == "single") {

          shinyjs::disable("l_param")

        } else {

          shinyjs::enable("l_param")

        }

      }

    }) # end controls setup


  })
}


#' @title Coloc package standard dataset
#' @description
#' Create the standard dataset for the [coloc] package
#' @param dat a data.table, with at least columns c("RSID","EA","OA","BP","BETA","SE") +/- c("EAF","N"),
#' @param type a string, either c("quant","cc")
#' @param sdY a numeric, std of the variable
#' @param N an integer, the sample size
#' @param S a numeric, proportion of cases
#' @param ld an LD matrix
#' @return a coloc dataset
#' @export
#'
make_coloc_dataset <- function(dat, type, sdY=NA_real_, N=NA_real_, S=NA_real_, ld=NULL) {

  # checks
  type <- match.arg(type, choices = c("quant","cc"))
  stopifnot("Standard names missing from the input to make a coloc dataset" = all(c("RSID","EA","OA","BP","BETA","SE") %in% names(dat)))
  if(type=="quant") {
    if(is.na(sdY) && !(all(c("EAF","N") %in% names(dat)) | all(c(!is.na(N), "EAF" %in% names(dat))))) {
      showNotification("`sdY` or data columns `EAF` and `N` (or input value `N`) must be provided", type="error")
      return(NULL)
    }
  }

  # if using SuSiE LD based method, use alleles for RSID
  if(!is.null(ld)) {

    if(!"LD_RSID_allele" %in% names(dat)) {
      showNotification("If running SuSiE finemapping please provide `dat` and `ld` output from genepi.utils::ld_matrix()", type="error")
      return(NULL)
    }

    # make multiallelic variants unique ids
    dat[, RSID := LD_RSID_allele]
  }

  # test for and deal with duplicate IDs
  if(sum(duplicated(dat$RSID), na.rm=TRUE) > 0) {
    dat <- dat[dat[, .I[which.min(P)], by=RSID]$V1]
    showNotification("Duplicate RSIDs found in source 1, taking the lowest P-value variant", type="warning")
  }

  # test for and deal with NAs
  if(sum(rowSums(is.na(dat[,list(RSID,BP,BETA,SE,EAF)]))) > 0) {

    na_rows <- rowSums(is.na(dat[,list(RSID,BP,BETA,SE,EAF)])) > 0
    dat <- dat[!na_rows, ]
    showNotification(paste0(sum(na_rows), " rows with NA values removed"), type="warning")

  }

  # create the base coloc object
  D1 <- list(
    snp      = dat$RSID,
    position = dat$BP,
    beta     = dat$BETA,
    varbeta  = dat$SE ^ 2,
    type     = type
  )

  # add additional parameters if provided
  if(all(c("EAF","N") %in% names(dat))) {

    D1$MAF <- dat$EAF
    D1$N   <- max(dat$N, na.rm=T)

  } else if("EAF" %in% names(dat) & !"N" %in% names(dat) & !is.na(N)) {

    D1$MAF <- dat$EAF
    D1$N   <- N

  } else {

    showNotification("`N` must be provided either as a parameter, or as a column in `dat`", type="error")
    return(NULL)

  }

  # proportion of sample that are cases
  if(!is.na(S)) {

    D1$s <- S

  }

  # add the st.dev of the Y variable
  if(!is.na(sdY)) {

    D1$sdY <- sdY

  }

  # add the LD matrix
  if(!is.null(ld)) {

    D1$LD <- ld[D1$snp, D1$snp]

  }

  # check the dataset
  tryCatch({
    if(!is.null(ld)) {
      coloc::check_dataset(D1, req="LD")
    } else {
      coloc::check_dataset(D1)
    }
  },
  error=function(e) {
    showNotification(paste0("Coloc dataset check failed - ", e), type="error", duration=10)
    return(NULL)
  })

  # return the object
  return(D1)
}






calc_credible_set <- function(results, dat_join_to, coverage, result_type="finemap.abf") {

  # checks
  result_type <- match.arg(result_type, choices = c("finemap.abf",
                                                    "finemap.signals",
                                                    "susie_rss"))

  # susie output
  if(result_type == "susie_rss") {

    result <- data.table::data.table(stack( susieR::susie_get_cs(res=results, coverage=coverage)$cs ))
    result[, group := sub("L([0-9]+)","\\1",ind)]
    result[, group := factor(group, levels=sort(unique(as.integer(group))))]
    result[, pip := results$pip[result$values]]
    result[, RSID_alleles := names(results$pip[result$values])]
    result[, RSID  := sub("_[ACTG]+_[ACTG]+$", "", RSID_alleles)]
    result[, index := pip == max(pip, na.rm=TRUE), by="group"]
    result <- result[, list(RSID,index,group)]

  # finemap output
  } else if(result_type == "finemap.abf") {

    result <- data.table::as.data.table(results)
    result[, RSID := snp]

    # order and get cumulative posterior probability
    result[order(SNP.PP, decreasing=TRUE), CUMSUM_PP := cumsum(SNP.PP)]

    # determine credible set - flag as coloc (i.e. when is the cumulative probability > 95% / whatever indicated)
    result[order(CUMSUM_PP), hit := ifelse(.I <= which(CUMSUM_PP >= coverage)[1], TRUE, FALSE)]

    # credible set factor
    result[, group := ifelse(hit, 1, NA)]
    result[, group := factor(group, levels=c(1))]
    result[, max_pp  := max(SNP.PP, na.rm=TRUE), by="group"]
    result[, group_n := .N, by="group"]
    result[, index   := data.table::fcase(is.na(group), NA,
                                          group_n == 1, TRUE,
                                          SNP.PP==max_pp,   TRUE,
                                          default = FALSE)]

  # finemap signals output
  } else if(result_type == "finemap.signals") {

    result <- data.table::data.table(RSID  = sub("_[ACTG]+_[ACTG]+$", "", names(results)),
                                     index = TRUE,
                                     group = factor(1:length(results)))

  #   # as data.table and clean up
  #   result <- data.table::as.data.table(results)
  #   data.table::setnames(result, c("SNP.PP") , c("SNP.PP1"))
  #   result <- result[snp != "null", ]
  #   result[, RSID := sub("_[ACTG]+_[ACTG]+$", "", snp)]
  #   # result[, RSID := snp]
  #
  #   # pivot longer
  #   result <- data.table::melt(result,
  #                              id.vars       = c("RSID"),
  #                              measure.vars  = grep("SNP.PP", names(result), value=TRUE),
  #                              variable.name = "STEP",
  #                              value.name    = "PP")
  #   result[, STEP := sub("SNP.PP([0-9]+)","\\1",STEP)]
  #
  #   # process groups
  #   result <- split(result, by="STEP")
  #
  #   result <- lapply(result, function(g) {
  #
  #     # order and get cumulative posterior probability
  #     g[order(PP, decreasing=TRUE), CUMSUM_PP := cumsum(PP)]
  #
  #     # determine credible set - flag as coloc (i.e. when is the cumulative probability > 95% / whatever indicated)
  #     g[order(CUMSUM_PP), hit := ifelse(.I <= which(CUMSUM_PP >= coverage)[1], TRUE, FALSE)]
  #
  #     # if not a assign group
  #     g[, group := factor(STEP, levels=sort(unique(as.integer(STEP))))]
  #     g[hit==FALSE, group := NA]
  #
  #   }) |> data.table::rbindlist()
  #
  #   # credible set factor
  #   result[, max_pp  := max(PP, na.rm=TRUE), by="group"]
  #   result[, group_n := .N, by="group"]
  #   result[, index   := data.table::fcase(is.na(group), NA,
  #                                         group_n == 1, TRUE,
  #                                         PP==max_pp,   TRUE,
  #                                         default = FALSE)]

    # browser()
    #
    # # only hits data
    # result <- result[hit==TRUE, ]
  }



  # join to the provided data
  d <- data.table::copy(dat_join_to) # force memory change to trigger reactives
  d[result, c("index","group") := list(i.index, i.group), on="RSID"]

  # return the data
  return(d)
}


#' @title Estimate trait variance, internal function
#' @description
#' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
#'
#' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
#' var(X) = 2*maf*(1-maf)
#' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
#'
#' @param vbeta vector of variance of coefficients
#' @param maf vector of MAF (same length as vbeta)
#' @param n sample size
#' @return estimated standard deviation of Y
#' @references https://github.com/chr1swallace/coloc/blob/a44c817137daa0f4d0de56c9b77693f8228ab5fc/R/claudia.R#L135C1-L157C2
#' @author Chris Wallace
#' @noRd
sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}


