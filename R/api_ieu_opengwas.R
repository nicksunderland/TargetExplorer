# API tutorial http://gwas-api.mrcieu.ac.uk/docs/
#' @title IEU Open GWAS API
#' @return a data.table of study information
#' @export
#'
get_ieu_gwas_info <- function() {

  # try to pull the dataset information from the API
  tryCatch({

    shiny::withProgress(message = 'Querying IEU OpenGWAS API', value = 0, {
      n=2
      shiny::incProgress(1/n, detail = "Getting study IDs")

      URL = "http://gwas-api.mrcieu.ac.uk/gwasinfo"

      # send the request
      r <- httr::GET(URL, httr::accept_json())
      httr::status_code(r)

      # get the content
      cont <- httr::content(r, "text", encoding = "UTF-8")

      # NULL if failed
      if (httr::status_code(r) != 200) {

        shiny::incProgress(1/n, detail = "Failed")
        showNotification(glue::glue("Error {httr::status_code(r)}"), type="error")
        return(NULL)

      }

      shiny::incProgress(1/n, detail = "Complete")
      gwas_list <- jsonlite::fromJSON(cont)
      df <- lapply(gwas_list, data.table::as.data.table) |> data.table::rbindlist(fill=TRUE)

      # return
      return(df)

    })

  }, error=function(e){

    print(e)
    return(NULL)

  })

}



#' @title Import IEU GWAS associations
#' @description
#' A function to download from the IEU OpenGWAS collection.
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
#'
import_ieu_gwas <- function(study_id,
                            chr      = NULL,
                            bp_start = NULL,
                            bp_end   = NULL,
                            build    = "GRCh37",
                            pthresh  = NULL) {

  # checks
  build <- match.arg(build, choices = c("GRCh37","GRCh38"))

  # lift position if needed
  if(build=="GRCh38") {
    pos_df <- data.frame(RSID=c("rs1","rs2"), CHR=c(chr,chr), BP=c(bp_start, bp_end))
    pos_df <- genepi.utils::lift(pos_df, from="Hg38", to="Hg19", snp_col="RSID", chr_col="CHR", pos_col="BP", remove_duplicates = FALSE)
    bp_start <- pos_df[1, ][["BP"]]
    bp_end   <- pos_df[2, ][["BP"]]
  }

  shiny::withProgress(message = 'Querying IEU OpenGWAS API', value = 0, {
    n=3

    # increment progress
    shiny::incProgress(1/n, detail = paste0("Dataset ", study_id, "... associations=", chr, ":", bp_start, "-", bp_end))

    # the base URL
    URL = glue::glue("http://gwas-api.mrcieu.ac.uk/associations")

    # send the request
    r <- httr::POST(url  = URL,
                    body = list(id      = study_id,
                                variant = paste0(chr,':',bp_start,'-',bp_end)),
                    httr::accept_json())
    cont <- httr::content(r, "text", encoding = "UTF-8")

    shiny::incProgress(1/n, detail = paste0("Dataset ", study_id, "... associations=", chr, ":", bp_start, "-", bp_end))

  }) # end withProgress


  # If the request was unsuccessful
  if (httr::status_code(r) != 200) {

    showNotification(glue::glue("Error {httr::status_code(r)}"), type="error")
    return(NULL)

  }

  # extract content
  data_out <- jsonlite::fromJSON(cont)

  # standardise
  if(build=="GRCh38") {
    data_out <- standardise_data(data_out, source="ieu_opengwas", build_from="Hg19", build_to="Hg38") # I think all of the GWAS in Open GWAS are b37
  } else if(build=="GRCh37") {
    data_out <- standardise_data(data_out, source="ieu_opengwas", build_from=NULL, build_to=NULL)
  }

  # apply pvalue threshold
  data_out <- data_out[P <= pthresh, ]

  # return the data
  return(data_out)

}
