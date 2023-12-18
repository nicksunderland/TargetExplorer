#' @title Ensembl gene information
#' @description
#' Query the ensembl REST API for gene information
#' @param gene_id a string, gene ID e.g. ENSG00000112164
#' @param build a string, one of c("GRCh37", "GRCh38")
#' @return a list of gene information - see https://grch37.rest.ensembl.org/lookup/id/ docs for structure
#' @export
#'
get_ensembl_gene_info <- function(gene_id, build) {

  # create the request URL
  if(build=="GRCh37") {

    url <- glue::glue("https://grch37.rest.ensembl.org/lookup/id/{trimws(gene_id)}?expand=1")

  } else if(build=="GRCh38") {

    url <- glue::glue("https://rest.ensembl.org/lookup/id/{trimws(gene_id)}?expand=1")

  }

  # send the request
  r <- httr::GET(url=url, config=httr::content_type("application/json"))

  # receive resonse and extract content
  httr::stop_for_status(r)
  gene_dat <- httr::content(r)

  return(gene_dat)

}
