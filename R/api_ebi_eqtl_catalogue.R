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
                                  study_info_table$quant_method=="ge", # TODO: not sure about this restriction/filter at the minute...
                                  c("dataset_id","catalogue","dataset")]
  dataset_ids <- unique(dataset_ids)

  # set up parameters - only done position and pvalue filtering for now
  parameters = list('pos'        = ifelse(all(!sapply(list(chr,bp_start,bp_end), is.null)), glue("{chr}:{bp_start}-{bp_end}"), ""),
                    'nlog10p'    = ifelse(nlog10p==0, "", nlog10p))

  # df to add results to
  data_out <- data.table::data.table()

  shiny::withProgress(message = 'Querying EBI eQTL API', value = 0, {

    # for progress bar
    n <- nrow(dataset_ids)*2*8 # allow roughly 8 hits per dataset

    # cycle the datasets
    for(id in dataset_ids$dataset_id) {
      # setup for this study
      size  = 1000
      start = 0

      # for each dataset, grab data in 1000 chuncks
      while (T) {
        URL = glue::glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{id}/associations?size={size}&start={start}")

        # increment progress
        shiny::incProgress(1/n, detail = paste0("Dataset ", id, "... ?size=", size, "&start=", start, "&pos=", parameters[["pos"]]))

        #Adding defined parameters to the request
        for (i in 1:length(parameters)) {
          if(is.null(parameters[[i]]) || parameters[[i]]=="") {
            URL = URL
          } else {
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
            print(glue::glue("Error {httr::status_code(r)}"))
            print(cont)
            return ()
          }
          #else just break
          break
        }

        # extract content
        cont_df <- jsonlite::fromJSON(cont)

        if (start == 0) {
          responses <- data.table::data.table(cont_df)
        } else {
          responses <- rbind(responses, cont_df)
        }
        start <- start + size
      }

      # add dataset info to the table
      responses[, STUDY  := dataset_ids[dataset_ids$dataset_id==id, "catalogue"]]
      responses[, TISSUE := dataset_ids[dataset_ids$dataset_id==id, "dataset"]]

      # add to the main table
      data_out <- rbind(data_out, responses)

      # update progress
      shiny::incProgress(1/n, detail = paste("Dataset", id, "complete"))

    }

  }) # end withProgress

  # No data found, return null
  if(nrow(data_out)==0) return(NULL)

  # calculate N as allele number / 2
  data_out[, N := ceiling(an/2)]

  # see description...
  if(!is.null(gene_id)) {
    data_out <- data_out[gene_id==gene_id, ]
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
    URL = glue::glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}&quant_method=ge")

    # send the request
    r <- httr::GET(URL, httr::accept_json())
    httr::status_code(r)

    # get the content
    cont <- httr::content(r, "text", encoding = "UTF-8")
    df <- jsonlite::fromJSON(cont)

    # standardise names
    names(df)[names(df)=="study_label"] <- "catalogue"
    names(df)[names(df)=="tissue_label"] <- "dataset"

    # return
    return(df)

  }, error=function(e){

    print(e)

  })

}
