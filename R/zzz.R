.onLoad <- function(libname, pkgname) {

  # set up caching for functions that hit the various APIs
  cache <- disk_cache()

  # study information to populate the search filters
  get_ebi_eqtl_info     <<- memoise::memoise(get_ebi_eqtl_info,     cache = cache)
  get_ebi_gwas_info     <<- memoise::memoise(get_ebi_gwas_info,     cache = cache)
  get_ieu_gwas_info     <<- memoise::memoise(get_ieu_gwas_info,     cache = cache)
  get_ensembl_gene_info <<- memoise::memoise(get_ensembl_gene_info, cache = cache)

  # data importing functions
  import_ieu_gwas       <<- memoise::memoise(import_ieu_gwas,       cache = cache)
  import_ebi_eqtl       <<- memoise::memoise(import_ebi_eqtl,       cache = cache)
  import_ebi_gwas       <<- memoise::memoise(import_ebi_gwas,       cache = cache)

}


#' @title Create a disk cache
#' @description
#' Returns a [cachem::cache_disk()] object with the cache file path located in this packages
#' `inst/cache` folder
#' @return a [cachem::cache_disk()] object
#' @importFrom cachem cache_disk
#' @export
#'
disk_cache <- function() {
  cache_path <- get_data_directory("cache")
  cache <- cachem::cache_disk(cache_path, max_size=10*1024*1024^2, max_age=60*60*24*90) # cache 10Gb and 90 day max
  return(cache)
}
