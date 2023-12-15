remove_shiny_inputs <- function(id, .input) {
  invisible(
    lapply(grep(id, names(.input), value=TRUE), function(i) {
      .subset2(.input, "impl")$.values$remove(i)
    })
  )
}


#' @title Get current operating system
#' @return a string, the oeprating system
#' @export
#'
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}




