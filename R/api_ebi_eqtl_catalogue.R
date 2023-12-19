#' @title Import EBI eQTL data
#' @description
#' A function to download from the EBI eQTL API
#' Note that there is some unexpected behaviour with the version2 filtering.
#' I thought that filtering by position and gene would be an `&` operation,
#' but it seems to be an `|` combination as when used together you get gene
#' eQTLs from way outside your specified region. Therefore I do the gene
#' filtering after.
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
#' @importFrom jsonlite fromJSON
#' @importFrom shiny withProgress incProgress
#'
import_ebi_eqtl <- function(study_info_table,
                            studies  = NULL,  #studies = c("BLUEPRINT","GTEx")
                            tissues  = NULL,  #tissues = c("adipose","adipose","monocyte")
                            method   = NULL,  #methods = c("ge", "exon" etc.)
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

  # lift position if needed
  if(build=="GRCh37") {
    pos_df <- data.frame(RSID=c("rs1","rs2"), CHR=c(chr,chr), BP=c(bp_start, bp_end))
    pos_df <- genepi.utils::lift(pos_df, from="Hg19", to="Hg38", snp_col="RSID", chr_col="CHR", pos_col="BP", remove_duplicates = FALSE)
    bp_start <- pos_df[1, ][["BP"]]
    bp_end   <- pos_df[2, ][["BP"]]
  }

  # work out which datasets we need to query
  dataset_ids <- study_info_table[study_info_table$catalogue  %in% studies &
                                  study_info_table$dataset %in% tissues &
                                  study_info_table$quant_method %in% method,
                                  c("dataset_id","catalogue","dataset","quant_method")]
  dataset_ids <- unique(dataset_ids)

  shiny::withProgress(message = 'Querying EBI eQTL API', value = 0, {

    # for progress bar
    n <- nrow(dataset_ids)

    # cycle the datasets
    data_out <- lapply(dataset_ids$dataset_id, function(id) {

      # extract the study/dataset (isolated function so that we can memoise cache the result)
      responses <- import_ebi_eqtl_dataset(id, chr, bp_start, bp_end)

      # no data found
      if(is.null(responses)) {

        shiny::incProgress(1/n, detail = paste("\nDataset", id, "failed"))
        return(NULL)

      # data found
      } else {

        # add dataset info to the table
        responses[, STUDY        := dataset_ids[dataset_ids$dataset_id==id, "catalogue"]]
        responses[, TISSUE       := dataset_ids[dataset_ids$dataset_id==id, "dataset"]]
        responses[, QUANT_METHOD := dataset_ids[dataset_ids$dataset_id==id, "quant_method"]]

        # update progress
        shiny::incProgress(1/n, detail = paste("\nDataset", id, "complete"))

        # return
        return(responses)
      }

    }) |> # end lapply datasets
      data.table::rbindlist()

  }) # end withProgress

  # No data found in any dataset, return null
  if(nrow(data_out)==0) return(NULL)

  # calculate N as allele number / 2
  data_out[, N := ceiling(an/2)]

  # see description...
  gene_id_input <- gene_id # same name as column, rename to avoid confusion in data.table
  if(!is.null(gene_id_input)) {
    data_out <- data_out[gene_id==gene_id_input, ]
  }

  # standardise
  if(build=="GRCh37") {
    data_out <- standardise_data(data_out, source="ebi_eqtl_api", build_from="Hg38", build_to="Hg19")
  } else if(build=="GRCh38") {
    data_out <- standardise_data(data_out, source="ebi_eqtl_api", build_from=NULL, build_to=NULL)
  }

  # return the data
  return(data_out)

}



#' @title Import EBI eQTL dataset
#' @param id the EBI eQTL study ID
#' @inheritParams import_ebi_eqtl
#' @return a data.table, the eQTL dataset with standardised column names
#' @export
#'
import_ebi_eqtl_dataset <- function(id, chr, bp_start, bp_end) {

  shiny::withProgress(message = paste0('Extracting dataset ', id), value = 0, {

    # setup for this study
    size  = 1000
    start = 0
    n = 8

    # for each dataset, grab data in 1000 chuncks
    while (T) {

      # increment progress
      shiny::incProgress(1/n, detail = glue::glue("size={size}&start={start}&pos={chr}:{bp_start}-{bp_end}"))

      # the URL - ENSG00000112164 - QTS000015 / QTD000141  &molecular_trait_id={gene_id} doesnt work as seems to cut the data off before the gene starts....
      URL = glue::glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{id}/associations?size={size}&start={start}&pos={chr}:{bp_start}-{bp_end}")

      # send the request
      r <- httr::GET(URL, httr::accept_json())
      cont <- httr::content(r, "text", encoding = "UTF-8")

      # If the request was unsuccessful
      if (httr::status_code(r) != 200) {
        if(start==0) {
          #If we get no results at all, print error
          showNotification(glue::glue("Error {httr::status_code(r)} - no results found in dataset={id}"), type="warning")
          return(NULL)
        } else {
          return(responses)
        }
      }

      # extract content
      cont_df <- jsonlite::fromJSON(cont) |> data.table::as.data.table()

      if (start == 0) {
        responses <- data.table::data.table(cont_df)
      } else {
        responses <- rbind(responses, cont_df)
      }

      # browser()

      start <- start + size
    }

  }) # end with progress

}






# API tutorial https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md
#' @title EBI eQTL Catalogue dataset information
#' @return a data.frame of study information
#' @export
#'
get_ebi_eqtl_info <- function() {

  # try to pull the dataset information from the API
  tryCatch({

    max_pulled_rows = 1000 # All datasets will be pulled if this parameter is bigger than the actual number of datasets

    # TODO: work out which quantification methods we want - gone with 'ge' gene counts, the default for now, but there are lots of different ones
    URL = glue::glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}") #&quant_method=ge

    # send the request
    r <- httr::GET(URL, httr::accept_json())
    httr::status_code(r)

    if (httr::status_code(r) != 200) {
      #If we get no results at all, print error
      showNotification(glue::glue("Status {httr::status_code(r)} - no EBI eQTL study information found"), type="error")
      return(NULL)
    }

    # get the content
    cont <- httr::content(r, "text", encoding = "UTF-8")
    df <- jsonlite::fromJSON(cont) |> data.table::as.data.table()

    # standardise names
    data.table::setnames(df, "study_label", "catalogue")
    data.table::setnames(df, "tissue_label", "dataset")

    # return
    return(df)

  }, error=function(e){

    print(e)
    return(NULL)

  })

}
