MR <- setClass(
  Class = "MR",
  slots = list(
    snps          = 'character',
    ea            = 'character',
    oa            = 'character',
    eafx          = 'matrix',
    nx            = 'matrix',
    bx            = 'matrix',
    bxse          = 'matrix',
    px            = 'matrix',
    by            = 'numeric',
    byse          = 'numeric',
    py            = 'numeric',
    eafy          = 'numeric',
    ny            = 'integer',
    exposure_id   = 'character',
    exposure      = 'character',
    outcome_id    = 'character',
    outcome       = 'character',
    correlation   = 'matrix'
  ),
  prototype = list(
    snps        = vector(mode='character'),
    ea          = vector(mode='character'),
    oa          = vector(mode='character'),
    eafx        = matrix(),
    nx          = matrix(),
    bx          = matrix(),
    bxse        = matrix(),
    px          = matrix(),
    by          = vector(mode='numeric'),
    byse        = vector(mode='numeric'),
    py          = vector(mode='numeric'),
    eafy        = vector(mode='numeric'),
    ny          = vector(mode='integer'),
    exposure    = vector(mode="character"),
    exposure_id = vector(mode='character'),
    outcome     = character(),
    outcome_id  = character(),
    correlation = matrix()
  )
)

setValidity("MR", function(object) {

  # correct lengths
  l <- list(object@snps, object@ea, object@oa, object@eafx[,1], object@bx[,1], object@bxse[,1], object@px[,1], object@by, object@byse, object@py)
  stopifnot("Inconsistent lengths" = all(lengths(l)==length(l[[1]])))

  # if LD matrix ensure symmetric and same size
  if(!all(is.na(object@correlation))) {
    stopifnot("LD matrix names do not match" = all(object@snps==rownames(object@correlation)))
    stopifnot("LD matrix doesn't match" = nrow(object@correlation)==length(l[[1]]))
  }

})


setGeneric("to_prune_correlation", function(ld_mat, prune_r2_thresh=0.95, seed=2024) standardGeneric("to_prune_correlation"))
setMethod("to_prune_correlation", "matrix", function(ld_mat, prune_r2_thresh=0.95, seed=2024) {

  # https://wellcomeopenresearch.org/articles/8-449
  # threshold set to r2>0.95
  set.seed(seed)        # for reproducibility
  omit = NULL           # set up list of variants to be omitted
  rho.upper = ld_mat    # correlation matrix
  rho.upper[lower.tri(ld_mat, diag=TRUE)] <- 0
  # only consider upper triangle of correlations
  j=1                   # set counter to 1

  while (max(abs(rho.upper), na.rm=TRUE) > prune_r2_thresh) {
    omit[j] = ifelse(rbinom(1, 1, 0.5)==1,
                     which.max(apply(abs(rho.upper), 1, max, na.rm=TRUE)),
                     which.max(apply(abs(rho.upper), 2, max, na.rm=TRUE)))
    # find the highest correlation value
    # select either the row or column at random
    # add this to the list of omitted variants
    rho.upper[omit[j],] <- 0 # set the correlations in this row to zero
    rho.upper[,omit[j]] <- 0 # set the correlations in this column to zero
    # (to avoid selecting the same variant again)
    j=j+1                # increment the counter
  }  # stop when no more pairwise correlations exceed threshold

  return(omit)
})

setGeneric("remove_variants", function(x, index, correlation_too=FALSE) standardGeneric("remove_variants"))
setMethod("remove_variants", "MR", function(x, index, correlation_too=FALSE) {

  # do the removing
  if(length(index)>0) {
    x@snps <- x@snps[-index]
    x@ea   <- x@ea[-index]
    x@oa   <- x@oa[-index]
    x@eafx <- as.matrix(x@eafx[-index, ])
    x@nx   <- as.matrix(x@nx[-index, ])
    x@bx   <- as.matrix(x@bx[-index, ])
    x@bxse <- as.matrix(x@bxse[-index, ])
    x@px   <- as.matrix(x@px[-index, ])
    x@by   <- x@by[-index]
    x@byse <- x@byse[-index]
    x@py   <- x@py[-index]
    x@eafy <- x@eafy[-index]
    x@ny   <- x@ny[-index]
    if(correlation_too) {
      x@correlation <- x@correlation[-index, -index]
    }
  }

  # validate and return
  validObject(x)
  return(x)
})

setGeneric("mr_load", function(object, exposure, outcome, e_map=NULL, o_map=NULL, harmonise_strictness=2) standardGeneric("mr_load"))
setMethod("mr_load", "MR", function(object, exposure, outcome, e_map=NULL, o_map=NULL, harmonise_strictness=2) {

  exp <- data.table::as.data.table(exposure)
  out <- data.table::as.data.table(outcome)

  #######
  #######
  ######
  # exp <- readRDS("/Users/xx20081/Desktop/exp.RDS")
  # out <- readRDS("/Users/xx20081/Desktop/out.RDS")
  # harmonise_strictness=2
  # e_map=NULL
  # o_map=NULL
  #######
  #######
  ######

  base_cols <- c('RSID','CHR','BP','EA','OA','EAF','BETA','SE','P','N','ID','TRAIT')
  names(base_cols) <- base_cols
  if(is.null(e_map)) { e_map <- base_cols }
  if(is.null(o_map)) { o_map <- base_cols }
  stopifnot("Missing mapping columns" = all(c(names(e_map), names(o_map)) %in% names(base_cols)))

  # map names
  data.table::setnames(exp, e_map, names(e_map))
  data.table::setnames(out, o_map, names(o_map))

  # IDs unique to allele coding
  exp[, sortedRSID := paste0(c(RSID,sort(c(OA,EA))), collapse="_"), by=1:nrow(exp)]
  out[, sortedRSID := paste0(c(RSID,sort(c(OA,EA))), collapse="_"), by=1:nrow(out)]

  # format function from TwoSampleMR
  format <- getFromNamespace('format_data', 'TwoSampleMR')

  # TwoSampleMR format args
  args <- list(snp_col           = "sortedRSID",
               chr_col           = "CHR",
               pos_col           = "BP",
               effect_allele_col = "EA",
               other_allele_col  = "OA",
               eaf_col           = "EAF",
               beta_col          = "BETA",
               se_col            = "SE",
               pval_col          = "P",
               samplesize_col    = "N",
               id_col            = "ID",
               phenotype_col     = "TRAIT")

  # format exposure - # rs2284795
  exp <- do.call(format, c(list(dat=as.data.frame(exp), type="exposure"), args))

  # format outcome
  out <- do.call(format, c(list(dat=as.data.frame(out), type="outcome"), args))

  # have the reference dataset as the first one in the table
  ref_dataset <- subset(exp, id.exposure == id.exposure[1])

  # TwoSampleMR::harmonise_data flips the outcome data, so have the full data as `outcome` type initially
  names(exp) <- gsub("exposure", "outcome", names(exp))

  # harmonise exposures against first exposure dataset
  harm_exp <- TwoSampleMR::harmonise_data(ref_dataset, exp, action=harmonise_strictness)
  harm_exp <- subset(harm_exp, mr_keep)

  # only keep harmonised SNPs common to all exposures
  tab <- table(harm_exp$SNP)
  keepsnps <- names(tab)[tab == length(unique(harm_exp$id.outcome))]
  harm_exp <- harm_exp[harm_exp$SNP %in% keepsnps, ]

  # drop dummy `exposure` columns and rename outcome to exposure
  harm_exp <- harm_exp[, !grepl("exposure|action|mr_keep|SNP_index", names(harm_exp))]
  names(harm_exp) <- gsub("outcome","exposure", names(harm_exp))

  # sanity checks before contoinuing and merging in allele data
  harm_exp <- data.table::as.data.table(harm_exp)
  harm_exp[, ea_check := all(effect_allele.exposure==effect_allele.exposure[[1]]), by="SNP"]
  harm_exp[, oa_check := all(other_allele.exposure==other_allele.exposure[[1]]), by="SNP"]
  stopifnot("Exposures harmonisation has failed" = all(c(harm_exp$ea_check, harm_exp$oa_check)))

  # harmonise the outcome to the same reference dataset
  harm_out <- TwoSampleMR::harmonise_data(harm_exp, out, action=harmonise_strictness)
  harm_out <- subset(harm_out, mr_keep)
  harm_out <- data.table::as.data.table(harm_out)

  # rename SNPs now harmonised
  harm_out[, SNP := sub("(.*)_(?:[ACTG]+|[DI])_(?:[ACTG]+|[DI])$","\\1",SNP,ignore.case=TRUE)]
  harm_out[, SNP := paste0(SNP,"_",other_allele.outcome,"_",effect_allele.outcome)]

  # beta matrix
  exp_beta  <- data.table::dcast(harm_out, SNP + effect_allele.exposure + other_allele.exposure ~ id.exposure, value.var='beta.exposure')
  object@bx <- as.matrix(exp_beta[,-c(1:3)])

  # SNP / alleles vectors
  object@snps <- exp_beta$SNP
  object@ea   <- exp_beta$effect_allele.exposure
  object@oa   <- exp_beta$other_allele.exposure

  # pvalue matrix
  exp_pval  <- data.table::dcast(harm_out, SNP ~ id.exposure, value.var='pval.exposure')
  object@px <- as.matrix(exp_pval[,-1])

  # sterr
  exp_se      <- data.table::dcast(harm_out, SNP ~ id.exposure, value.var='se.exposure')
  object@bxse <- as.matrix(exp_se[,-1])

  # eaf
  exp_eaf     <- data.table::dcast(harm_out, SNP ~ id.exposure, value.var='eaf.exposure')
  object@eafx <- as.matrix(exp_eaf[,-1])

  # sample size
  exp_n     <- data.table::dcast(harm_out, SNP ~ id.exposure, value.var='samplesize.exposure')
  object@nx <- as.matrix(exp_n[,-1])

  # outcome vectors
  object@by   <- harm_out[id.exposure==id.exposure[[1]], beta.outcome]
  object@byse <- harm_out[id.exposure==id.exposure[[1]], se.outcome]
  object@py   <- harm_out[id.exposure==id.exposure[[1]], pval.outcome]
  object@eafy <- harm_out[id.exposure==id.exposure[[1]], eaf.outcome]
  object@ny   <- as.integer(harm_out[id.exposure==id.exposure[[1]], samplesize.outcome])

  # trait and id
  object@exposure    <- harm_out[, unique(exposure), by=id.exposure]$V1
  object@exposure_id <- harm_out[, unique(id.exposure), by=id.exposure]$V1
  object@outcome     <- harm_out$outcome[[1]]
  object@outcome_id  <- harm_out$id.outcome[[1]]

  # validate and return
  validObject(object)
  return(object)
})

setGeneric("is_multivariable", function(x) standardGeneric("is_multivariable"))
setMethod("is_multivariable", "MR", function(x) { ncol(x@bx)>1 })
setGeneric("num_exposures", function(x) standardGeneric("num_exposures"))
setMethod("num_exposures", "MR", function(x) { ncol(x@bx) })
setGeneric("snps_no_alleles", function(x) standardGeneric("snps_no_alleles"))
setMethod("snps_no_alleles", "MR", function(x) {
  return( sub("(.*)_(?:[ACTG]+|[DI])_(?:[ACTG]+|[DI])$", "\\1", x@snps, ignore.case=TRUE) )
})
setGeneric("set_ld_mat", function(x, ld_mat, rsid, ref, alt, align_ref_ea=TRUE, prune_r2_thresh=0.95) standardGeneric("set_ld_mat"))
setMethod("set_ld_mat", "MR", function(x, ld_mat, rsid, ref, alt, align_ref_ea=TRUE, prune_r2_thresh=0.95) {

  # ########
  # ld_mat = readRDS("/Users/xx20081/Desktop/ld_mat.RDS")
  # rsid = readRDS("/Users/xx20081/Desktop/rsid.RDS")
  # ref = readRDS("/Users/xx20081/Desktop/ref.RDS")
  # alt = readRDS("/Users/xx20081/Desktop/alt.RDS")
  # x=readRDS("/Users/xx20081/Desktop/x.RDS")
  #########

  # current (harmonised) alleles
  current <- data.table::data.table(RSID=snps_no_alleles(x), EA=x@ea, OA=x@oa)
  ld_alleles <- data.table::data.table(RSID=rsid, REF=ref, ALT=alt)[, RSID_alleles := paste0(RSID,"_",REF,"_",ALT)]

  # join and flag those that need flipping
  current[ld_alleles, LD_flip := !align_ref_ea, on=c("RSID"="RSID","EA"="REF","OA"="ALT")]
  current[ld_alleles, LD_flip :=  align_ref_ea, on=c("RSID"="RSID","EA"="ALT","OA"="REF")]

  # the indices of those that need flipping and those that need removing
  flip <- which(current$LD_flip==TRUE)
  nold <- which(is.na(current$LD_flip))

  # do the flipping
  ea_store       <- x@ea
  x@ea[flip]     <- x@oa[flip]
  x@oa[flip]     <- ea_store[flip]
  x@snps         <- paste0(snps_no_alleles(x),"_",x@ea,"_",x@oa)
  x@eafx[flip, ] <- as.matrix(1-x@eafx[flip, ])
  x@bx[flip, ]   <- as.matrix(x@bx[flip, ]*-1)
  x@by[flip]     <- x@by[flip]*-1
  x@eafy[flip]   <- 1-x@eafy[flip]

  # do the removing
  x <- remove_variants(x, index=nold, correlation_too=FALSE)

  # ensure order of ld matrix matches data
  ld_mat <- ld_mat[x@snps, x@snps]

  # set the LD matrix
  x@correlation <- ld_mat

  # prune the matrix and remove pruned variants from object and matrix
  prune_idx <- to_prune_correlation(ld_mat, prune_r2_thresh=prune_r2_thresh)
  x <- remove_variants(x, index=prune_idx, correlation_too=TRUE)

  # validate and return
  validObject(x)
  return(x)
})


setGeneric("mr_results_to_plotting", function(x, id_idx=1, orientate=TRUE) standardGeneric("mr_results_to_plotting"))
setMethod("mr_results_to_plotting", "MR", function(x, id_idx=1, orientate=TRUE) {
  points_df <- data.frame(
    snp = x@snps,
    x   = x@bx[,id_idx],
    y   = x@by,
    xse = x@bxse[,id_idx],
    yse = x@byse
  )
  if(orientate) {
    neg_x <- points_df$x < 0
    points_df$x[neg_x] <- points_df$x[neg_x] * -1
    points_df$y[neg_x] <- points_df$y[neg_x] * -1
  }
  return(points_df)
})
setMethod("mr_results_to_plotting", "list", function(x) {

  mr_results <- lapply(x, function(r) parse_mr_result(r)) |> data.table::rbindlist(idcol="Method")

})


setGeneric("parse_mr_result", function(x) standardGeneric("parse_mr_result"))
setMethod("parse_mr_result", "IVW", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                SNPs     = x@SNPs,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue,
                Fstat    = x@Fstat)
})
setMethod("parse_mr_result", "MVIVW", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                SNPs     = x@SNPs,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue,
                CondFstat= x@CondFstat)
})
setMethod("parse_mr_result", "Egger", function(x) {
    mr_result_row(Exposure = x@Exposure,
                  Outcome  = x@Outcome,
                  SNPs     = x@SNPs,
                  Estimate = x@Estimate,
                  StdError = x@StdError.Est,
                  Pvalue   = x@Pvalue.Est,
                  Intercept= x@Intercept,
                  Int.SE   = x@StdError.Int,
                  Int.Pval = x@Pvalue.Int)
})
setMethod("parse_mr_result", "MVEgger", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                SNPs     = x@SNPs,
                Estimate = x@Estimate,
                StdError = x@StdError.Est,
                Pvalue   = x@Pvalue.Est,
                Intercept= x@Intercept,
                Int.SE   = x@StdError.Int,
                Int.Pval = x@Pvalue.Int)
})
setMethod("parse_mr_result", "WeightedMedian", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                SNPs     = x@SNPs,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue)
})
setMethod("parse_mr_result", "MVMedian", function(x) {

  print(x)

  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                SNPs     = x@SNPs,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue)
})
setMethod("parse_mr_result", "MRMBE", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                SNPs     = x@SNPs,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue)
})
setMethod("parse_mr_result", "PCGMM", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue,
                Fstat    = x@Fstat,
                Overdispersion = x@Overdispersion,
                PCs      = x@PCs)
})
setMethod("parse_mr_result", "MVPCGMM", function(x) {
  mr_result_row(Exposure = x@Exposure,
                Outcome  = x@Outcome,
                Estimate = x@Estimate,
                StdError = x@StdError,
                Pvalue   = x@Pvalue,
                CondFstat= x@CondFstat,
                Overdispersion = x@Overdispersion,
                PCs      = x@PCs)
})
mr_result_row <- function(Exposure = NA_character_,
                          Outcome  = NA_character_,
                          SNPs     = NA_integer_,
                          Estimate = NA_real_,
                          StdError = NA_real_,
                          Pvalue   = NA_real_,
                          Intercept= 0,
                          Int.SE   = NA_real_,
                          Int.Pval = NA_real_,
                          Fstat    = NA_real_,
                          CondFstat= NA_real_,
                          Overdispersion = NA_real_,
                          PCs      = NA_integer_) {

  s <- list(
      Exposure = as.character(Exposure),
      Outcome  = as.character(Outcome),
      SNPs     = as.integer(SNPs),
      Estimate = as.numeric(Estimate),
      StdError = as.numeric(StdError),
      Pvalue   = as.numeric(Pvalue),
      Intercept= as.numeric(Intercept),
      Int.SE   = as.numeric(Int.SE),
      Int.Pval = as.numeric(Int.Pval),
      Fstat    = as.numeric(Fstat),
      CondFstat= as.numeric(CondFstat),
      Overdispersion = as.numeric(Overdispersion),
      PCs      = as.integer(PCs)
    )
  max_len <- max(lengths(s))
  s <- lapply(s, function(x) {

    if(length(x)<max_len) {
      x <- rep(x[[1]], max_len)
    } else {
      x
    }
  })
  s$id_idx <- 1:max_len
  return(data.frame(do.call(cbind.data.frame, s)))
}


setGeneric("to_MRInput", function(x) standardGeneric("to_MRInput"))
setMethod("to_MRInput", "MR", function(x) {

  mr_input <- MendelianRandomization::mr_input(
    bx          = x@bx[,1],
    bxse        = x@bxse[,1],
    by          = x@by,
    byse        = x@byse,
    exposure    = x@exposure[1],
    outcome     = x@outcome,
    snps        = x@snps,
    correlation = x@correlation
  )

  return(mr_input)
})

setGeneric("to_MRMVInput", function(x) standardGeneric("to_MRMVInput"))
setMethod("to_MRMVInput", "MR", function(x) {

  mr_input <- MendelianRandomization::mr_mvinput(
    bx          = x@bx,
    bxse        = x@bxse,
    by          = x@by,
    byse        = x@byse,
    exposure    = x@exposure,
    outcome     = x@outcome,
    snps        = x@snps,
    correlation = x@correlation
  )

  return(mr_input)
})



setGeneric("mr_ivw", function(x) standardGeneric("mr_ivw"))
setMethod("mr_ivw", "MR", function(x) {
  if(is_multivariable(x)) {
    res <- MendelianRandomization::mr_mvivw(to_MRMVInput(x), nx=apply(x@nx,2,max,na.rm=TRUE))
  } else {
    res <- MendelianRandomization::mr_ivw(to_MRInput(x))
  }
  return(res)
})

setGeneric("mr_egger", function(x) standardGeneric("mr_egger"))
setMethod("mr_egger", "MR", function(x) {
  if(is_multivariable(x)) {
    res <- MendelianRandomization::mr_mvegger(to_MRMVInput(x))
  } else {
    res <- MendelianRandomization::mr_egger(to_MRInput(x))
  }
  return(res)
})

setGeneric("mr_weighted_median", function(x) standardGeneric("mr_weighted_median"))
setMethod("mr_weighted_median", "MR", function(x) {
  if(is_multivariable(x)) {
    res <- MendelianRandomization::mr_mvmedian(to_MRMVInput(x))
  } else {
    res <- MendelianRandomization::mr_median(to_MRInput(x), weighting="weighted")
  }
  return(res)
})

setGeneric("mr_weighted_mode", function(x) standardGeneric("mr_weighted_mode"))
setMethod("mr_weighted_mode", "MR", function(x) {
  if(is_multivariable(x)) {
    warning("No multivariable mode based function")
    return(NULL)
  } else {
    res <- MendelianRandomization::mr_mbe(to_MRInput(x), weighting='weighted')
  }
  return(res)
})

setGeneric("mr_pcgmm", function(x) standardGeneric("mr_pcgmm"))
setMethod("mr_pcgmm", "MR", function(x) {
  if(is_multivariable(x)) {
    res <- MendelianRandomization::mr_mvpcgmm(to_MRMVInput(x), nx=apply(x@nx,2,max,na.rm=TRUE), ny=max(x@ny,na.rm=TRUE))
  } else {
    res <- MendelianRandomization::mr_pcgmm(to_MRInput(x), nx=max(x@nx[,1],na.rm=TRUE), ny=max(x@ny,na.rm=TRUE))
  }
  return(res)
})






# # results_uv = readRDS("/Users/xx20081/Desktop/uni_results.RDS")
# # results_mv = readRDS("/Users/xx20081/Desktop/mv_results.RDS")
#
# exp    <- readRDS("/Users/xx20081/Desktop/exp.RDS")
# out    <- readRDS("/Users/xx20081/Desktop/out.RDS")
# object <- MR()
# object <- mr_load(object, exposure=exp, outcome=out)
#
#
# library(ggplot2)
# ggplot(data = mr_results_to_plotting(object), mapping = aes(x=x, y=y)) +
#   geom_errorbar( mapping = aes(ymin=y-yse, ymax=y+yse), width=0, color="grey") +
#   geom_errorbarh(mapping = aes(xmin=x-xse, xmax=x+xse), height=0, color="grey") +
#   geom_point()


# is_multivariable(object)
# num_exposures(object)
# res_ivw <- mr_ivw(object)
# res_egg <- mr_egger(object)
# # res_med <- mr_weighted_median(object)
# res_mod <- mr_weighted_mode(object)

