#' @title Import EBI GWAS associations
#' @description
#' A function to download from the EBI GWAS Catalogue. Build GRCh38 dbSNP Build 154
#' @param study_id .
#' @param chr .
#' @param bp_start .
#' @param bp_end .
#' @param build .
#' @param pthresh .
#' @return a data.frame
#' @export
#' @importFrom glue glue
#' @importFrom jsonlite fromJSON
#' @importFrom shiny withProgress incProgress
#' @importFrom cli cli_alert_info
#'
import_ebi_gwas <- function(study_id,
                            chr      = NULL,
                            bp_start = NULL,
                            bp_end   = NULL,
                            build    = "GRCh38",
                            pthresh  = NULL) {

  # checks
  build <- match.arg(build, choices = c("GRCh37","GRCh38"))

  # lift position if needed
  if(build=="GRCh37") {
    pos_df <- data.frame(RSID=c("rs1","rs2"), CHR=c(chr,chr), BP=c(bp_start, bp_end))
    pos_df <- genepi.utils::lift(pos_df, from="Hg19", to="Hg38", snp_col="RSID", chr_col="CHR", pos_col="BP", remove_duplicates = FALSE)
    bp_start <- pos_df[1, ][["BP"]]
    bp_end   <- pos_df[2, ][["BP"]]
  }

  # standard information format; contain to collect to
  data_template <- data.table::data.table(variant_id              = character(),
                                          chromosome              = character(),
                                          base_pair_location      = integer(),
                                          effect_allele           = character(),
                                          other_allele            = character(),
                                          effect_allele_frequency = character(),
                                          beta                    = numeric(),
                                          se                      = numeric(),
                                          odds_ratio              = numeric(),
                                          ci_lower                = numeric(),
                                          ci_upper                = numeric(),
                                          p_value                 = numeric(),
                                          trait                   = character())

  shiny::withProgress(message = 'Querying EBI GWAS Catalogue', value = 0, {
    n=8

    # increment progress
    shiny::incProgress(1/n, detail = paste0("Dataset ", study_id, "... associations=", chr, ":", bp_start, "-", bp_end))

    # the base URL
    URL = glue::glue("https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chr}/associations?")
    URL = glue::glue("{URL}study_accession={study_id}&size=1000&bp_lower={bp_start}&bp_upper={bp_end}&p_lower=0&p_upper=1")

    data_out <- data.table::data.table()
    start = 0

    # for each dataset, grab data in 1000 chuncks
    while (!is.null(URL)) {
      # increment progress
      shiny::incProgress(1/n)

      # send the request
      r <- httr::GET(URL, httr::accept_json())
      cont <- httr::content(r, "text", encoding = "UTF-8")

      # If the request was unsuccessful
      if (httr::status_code(r) != 200) {
        #If we get no results at all, print error
        if (start == 0) {
          print(glue::glue("Error {httr::status_code(r)}"))
          print(cont)
          return (NULL)
        }
        #else just break
        break
      }

      # extract content to a list of data.table rows
      cont_lst <- jsonlite::fromJSON(cont)
      cont_df <- lapply(cont_lst$`_embedded`$associations, function(row){
        row$`_links` <- NULL
        row <- data.table::as.data.table(row)
        row <- data.table::rbindlist(list(data_template, row), fill=TRUE)
        return(row)
      }) |>
        # bind rows into a table, based on the template columns; fill NAs
        data.table::rbindlist()

      # update URL to next page of results
      URL <- cont_lst$`_links`$`next`$href
      if(!is.null(URL)) {
        cli::cli_alert_info(paste0("Next URL: ", URL))
      }

      # flag moving to next page
      start <- start + 1

      # bind in to overall output
      data_out <- rbind(data_out, cont_df)

      # update progress
      shiny::incProgress(1/n, detail = paste0("Dataset ", study_id, "... associations=", chr, ":", bp_start, "-", bp_end))

    }

    # update progress
    shiny::incProgress(1/n, detail = paste("Dataset", study_id, "complete"))

  }) # end withProgress

  # standardise
  if(build=="GRCh37") {
    data_out <- standardise_data(data_out, source="ebi_gwas", build_from="Hg38", build_to="Hg19") # EBI is build 38
  } else if(build=="GRCh38") {
    data_out <- standardise_data(data_out, source="ebi_gwas", build_from=NULL, build_to=NULL)
  }

  # apply pvalue threshold
  data_out <- data_out[P <= pthresh, ]

  # return the data
  return(data_out)

}


# API docs https://www.ebi.ac.uk/gwas/rest/docs/api
#' @title EBI GWAS Catalogue dataset information
#' @return a data.frame of study information
#' @export
#'
get_ebi_gwas_info <- function() {

  withProgress(message = "Updating EBI GWAS Catalogue...", {
    n=3
    shiny::incProgress(1/n, detail = "contacting API")

    # base URL - request 100,000 (way more than in db)
    url <- "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies?size=1000" #  "https://www.ebi.ac.uk/gwas/rest/api/studies?page=0&size=500"

    # retrieve the page
    r <- httr::GET(url, httr::accept_json())
    httr::status_code(r)

    # If the request was unsuccessful
    if (httr::status_code(r) != 200) {
      #If we get no results at all, print error
      print(glue::glue("Error {httr::status_code(r)}"))
      print(cont)
      return (NULL)
    }

    cont <- httr::content(r, "text", encoding = "UTF-8")
    res  <- jsonlite::fromJSON(cont)

    # to table, with column name `id`
    shiny::incProgress(1/n, detail = "processing data")
    info <- data.table::data.table(id = sapply(res$`_embedded`$studies, function(x){ x$study_accession }))

    # return
    shiny::incProgress(1/n, detail = "complete")
    return(info)

  })

}
