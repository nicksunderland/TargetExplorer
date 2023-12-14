#' @title Import EBI eQTL data
#' @param study_info_table .
#' @param studies .
#' @param tissues .
#' @param chr .
#' @param bp_start .
#' @param bp_end .
#' @param cptid .
#' @param rsid .
#' @param gene_id .
#' @param nlog10p .
#' @param build .
#' @return a data.frame
#' @export
#' @importFrom glue glue
#' @importFrom shiny withProgress incProgress
#'
import_ebi_eqtl <- function(study_info_table,
                            studies  = NULL,  #studies = c("BLUEPRINT","GTEx")
                            tissues  = NULL,  #tissues = c("adipose","adipose","monocyte")
                            chr      = NULL,
                            bp_start = NULL,
                            bp_end   = NULL,
                            cptid    = NULL,
                            rsid     = NULL,
                            gene_id  = NULL,
                            nlog10p  = NULL,
                            build    = "GRCh37") {

  # checks
  build <- match.arg(build, choices = c("GRCh37","GRCh38"))

  # work out which datasets we need to query
  dataset_ids <- study_info_table[study_info_table$study_label  %in% studies &
                                  study_info_table$tissue_label %in% tissues &
                                  study_info_table$quant_method=="ge", # TODO: not sure about this restriction/filter at the minute...
                                  c("dataset_id","study_label","tissue_label")]
  dataset_ids <- unique(dataset_ids)

  # set up parameters

  parameters = list('pos'        = ifelse(all(!sapply(list(chr,bp_start,bp_end), is.null)), glue("{chr}:{bp_start}-{bp_end}"), ""),
                    'variant'    = cptid,
                    'rsid'       = rsid,
                    'gene_id'    = gene_id,
                    'nlog10p'    = nlog10p)

  # df to add results to
  data_out <- data.frame()

  shiny::withProgress(message = 'Querying EBI eQTL API', value = 0, {
    n <- nrow(dataset_ids)*2

    # cycle the datasets
    for(id in dataset_ids$dataset_id) {
      shiny::incProgress(1/n, detail = paste("Dataset", id, "..."))

      # setup for this study
      size  = 1000
      start = 0

      # for each dataset, grab data in 1000 chuncks
      while (T) {
        URL = glue::glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{id}/associations?size={size}&start={start}")

        #Adding defined parameters to the request
        for (i in 1:length(parameters)) {
          if(!is.null(parameters[[i]])) {
            URL = glue::glue("{URL}&{names(parameters[i])}={parameters[[i]]}")
          }
        }

        # send the request
        r <- httr::GET(URL, httr::accept_json())
        cont <- httr::content(r, "text", encoding = "UTF-8")

        # If the request was unsuccessful
        if (httr::status_code(r) != 200) {
          #If we get no results at all, print error
          if (start == 0) {
            print(glue::glue("Error {status_code(r)}"))
            print(cont)
            return ()
          }
          #else just break
          break
        }

        # extract content
        cont_df <- jsonlite::fromJSON(cont)

        if (start == 0) {
          responses <- cont_df
        } else {
          responses <- rbind(responses, cont_df)
        }
        start <- start + size
      }

      # add dataset info to the table
      responses$study_label  <- dataset_ids[dataset_ids$dataset_id==id, "study_label"]
      responses$tissue_label <- dataset_ids[dataset_ids$dataset_id==id, "tissue_label"]

      # add to the main table
      data_out <- rbind(data_out, responses)

      # update progress
      shiny::incProgress(1/n, detail = paste("Dataset", id, "complete"))

    }

  }) # end withProgress



  # standardise naming



  # lift to build 37 if needed as EBI catalogue is build 38
  if(build=="GRCh37") {
    data_out <- genepi.utils::lift(data_out, from="Hg38", to="Hg19",
                                   snp_col="rsid", chr_col="chromosome", pos_col="position", ea_col="alt", oa_col="ref",
                                   remove_duplicates = FALSE)
  }

  # return the data
  return(data_out)

}






















# API tutorial https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md
#' @title QBI eQTL Catalogue dataset information
#' @return a data.frame of study information
#' @export
#'
get_ebi_study_info <- function() {

  # try to pull the dataset information from the API
  tryCatch({

    max_pulled_rows = 1000 # All datasets will be pulled if this parameter is bigger than the actual number of datasets

    # TODO: work out which quantification methods we want - gone with 'ge' gene counts, the default for now, but there are lots of different ones
    URL = glue::glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}&quant_method=ge")

    # send the request
    r <- httr::GET(URL, httr::accept_json())
    httr::status_code(r)

    # get the content
    cont <- httr::content(r, "text", encoding = "UTF-8")
    df <- jsonlite::fromJSON(cont)

    # return
    return(df)

  }, error=function(e){

    print(e)

  })

}
