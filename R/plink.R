get_plink2_exe <- function() {

  # the plink2 executable
  if(get_os() == "osx") {
    plink2 <- system.file("resources", "plink2", package="TargetExplorer")
  } else if(get_os() == "linux") {
    plink2 <- system.file("resources", "plink2_linux", package="TargetExplorer")
  } else {
    message("Operating system: ", get_os())
    stop("operating system error")
  }

  # give the plink2 executable permission to run
  Sys.chmod(plink2, "777", use_umask = FALSE)

  return(plink2)
}
