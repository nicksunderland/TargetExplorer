make_1000G_ref_subset <- function(chrom, from, to) {

  # create file name for the reference subset
  cache_path <- file.path(get_data_directory("references"), "plink_1000gp3", "cache")
  out <- file.path(cache_path, paste0("ref1000gp3_chr",chrom,"_",from,"_",to))

  # see if it already exists, if so return it
  if(paste0(out,".pgen") %in% list.files(cache_path, full.names=TRUE)) {
    return(out)
  }

  # else, create a subset of the reference file
  plink2 <- get_plink2_exe()
  plink2_ref <- file.path(get_data_directory("references"), "plink_1000gp3", "all_phase3_nodup.pgen")

  # build the command
  cmd <- paste(plink2,
               "--pfile", sub(".pgen", "", plink2_ref),
               "--chr", chrom,
               "--from-kb", floor(from/1000),
               "--to-kb", ceiling(to/1000),
               "--make-pgen",
               "--out", out)

  # execute
  system(cmd)

  # confirm file
  if(!file.exists(paste0(out,".pgen"))) {
    stop("problems creating plink2 reference subset")
  } else {
    return(out)
  }
}


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
