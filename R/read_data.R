get_data_directory <- function(type=NULL) {

  # First look in the config file
  config <- yaml::read_yaml(system.file("app-config.yml", package="TargetExplorer"))
  if(dir.exists(config[["data_directory"]])) {
    if(is.null(type)) return(data_dir) else return(file.path(data_dir,type))
  }

  # look for installed data directory
  data_dir <- system.file("extdata", package="TargetExplorerData")
  if(dir.exists(data_dir)) {
    if(is.null(type)) return(data_dir) else return(file.path(data_dir,type))
  }

  # Look next to this packages directory
  this_package_dir <- dirname( dirname( system.file("", package="TargetExplorer") ) )
  data_dir <- file.path(this_package_dir, "TargetExplorerData", "inst", "extdata")
  if(dir.exists(data_dir)) {
    if(is.null(type)) return(data_dir) else return(file.path(data_dir,type))
  }

  stop("Data directory not found - have you installed TargetExplorerData?")
}


get_available_data_sources <- function(type) {

  if(type=="gwas"){

    # the each source is a directory containing the gwas file split into chromosomes (e.g. ".../hba1c_jurgens2022/chr1.RDS")
    sources <- list.files(get_data_directory(type), full.names=TRUE, recursive=FALSE)
    names(sources) <- list.files(get_data_directory(type), full.names=FALSE, recursive=FALSE)
    return(sources)

  } else if(type=="eqtl") {
    # TODO
  } else if(type=="pqtl") {
    # TODO
  }


}


read_internal_data <- function(type, source, chrom, start, end) {

  # data directory
  data_dir <- get_data_directory(type)

  # the data source files
  source_dir <- file.path(data_dir, source)

  # file path
  source_file <- file.path(source_dir, paste0("chr",chrom,".fst"))

  # create a fst object which allows row access without reading the whole file
  chr_fst <- fst::fst(source_file)

  # read just the position data (an integer vector)
  pos <- chr_fst[, "BP"]

  # apply the BP filter and get the indices
  row_idxs <- which(pos >= start & pos <= end)

  # read the needed rows
  d <- chr_fst[row_idxs, ]

  return(d)
}


read_gene_data <- function(source, chrom, start, end) {

  # data directory
  data_dir <- get_data_directory("genes")

  # the data source files
  source_dir <- file.path(data_dir, source)

  # file path
  source_file <- file.path(source_dir, paste0("chr",chrom,".fst"))

  # create a fst object which allows row access without reading the whole file
  chr_fst <- fst::fst(source_file)

  # read just the position data (an integer vector)
  pos <- chr_fst[, c("BP_START","BP_END")]

  # apply the BP filter and get the indices
  row_idxs <- which(pos$BP_END >= start & pos$BP_START <= end)

  # read the needed rows
  d <- chr_fst[row_idxs, ]

  # only one row per gene, for now
  d <- d[!duplicated(d$GENE_NAME), ]

  return(d)
}
