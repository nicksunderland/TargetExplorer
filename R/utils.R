remove_shiny_inputs <- function(id, .input, .session) {

  # remove all the inputs with a matching id
  module_inputs <- grep(id, names(.input), value=TRUE)
  lapply(module_inputs, function(i) { .subset2(.input, "impl")$.values$remove(i) })

  # destroy() then remove all the observers with a matching id
  module_observers <- grep(id, names(.session$userData), value=TRUE)
  lapply(module_observers, function(i) { .session$userData[[i]]$destroy() })
  rm(list=ls(module_observers, envir=.session$userData), envir=.session$userData)

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





