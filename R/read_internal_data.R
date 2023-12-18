get_data_directory <- function(type=NULL) {

  tryCatch({
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

  },error=function(e){
    print(e)
    # print(getwd())
    # return("/Users/xx20081/git/TargetExplorerData/extdata") # hack to pass devtools::check.......
  })


}


get_available_data_sources <- function() {

  source_info <- data.frame(catalogue=character(), dataset=character(), path=character())

  for(type in c("gwas","eqtl","pqtl","genes")) {

    s <- data.frame(path      = list.files(get_data_directory(type), full.names=TRUE,  recursive=FALSE),
                    dataset   = list.files(get_data_directory(type), full.names=FALSE, recursive=FALSE))

    if(nrow(s)>0) {
      s$catalogue <- type
      source_info <- rbind(source_info, s)
    }

  }

  return(source_info)
}


read_internal_data <- function(type, dataset, chrom, start, end) {

  shiny::withProgress(message = 'Reading package data', value = 0, {
    shiny::incProgress(1/4, detail = paste("Dataset", dataset, "..."))

    # data directory
    data_dir <- get_data_directory(type)

    # the data source files
    source_dir <- file.path(data_dir, dataset)

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

    # standardise data.frame
    d <- standardise_data(d, source="internal")

    shiny::incProgress(3/4, detail = paste("Complete"))
  })

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
