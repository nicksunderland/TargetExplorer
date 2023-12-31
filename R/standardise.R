#' @title Standardise data inputs for processing
#' @description
#' Minimum dataset for GWAS: c(RSID, CHR, BP, EA, OA, EAF, BETA, SE, P, N)
#' @param dat a data.frame like object
#' @param source a string, one of c("ebi_eqtl_api")
#' @param build_from valid parameter for `genepi.utils::lift()`
#' @param build_to valid parameter for `genepi.utils::lift()`
#' @return a data.table
#' @export
#'
standardise_data <- function(dat, source, build_from=NULL, build_to=NULL) {

  # R CMD checks
  nlog10P <- P <- NULL

  # checks
  source <- match.arg(source, choices = c("ebi_eqtl_api","internal","ieu_opengwas","ebi_gwas"))
  if(is.null(dat)) return(NULL)

  # work with data.tables
  dat <- data.table::as.data.table(dat)

  # EBI eQTL API - dealing with output from `import_ebi_eqtl()`
  if(source=="ebi_eqtl_api") {

    col_map <- c(RSID="rsid", CHR="chromosome", BP="position", EA="alt", OA="ref", EAF="maf", BETA="beta",SE="se", P="pvalue", N="N",
                 GENE_ID="gene_id", STUDY="STUDY", TISSUE="TISSUE", QUANT_METHOD="QUANT_METHOD", TPM="median_tpm",TRAIT="TRAIT", ID="ID")

  } else if(source=="internal") {

    col_map <- c(RSID="RSID", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", BETA="BETA", SE="SE", P="P", N="N",TRAIT="TRAIT", ID="ID")

  } else if(source=="ieu_opengwas") {

    col_map <- c(RSID="rsid", CHR="chr", BP="position", EA="ea", OA="nea", EAF="eaf", BETA="beta", SE="se", P="p", N="n", TRAIT="trait", ID="id")

  } else if(source=="ebi_gwas") {

    col_map <- c(RSID="variant_id", CHR="chromosome", BP="base_pair_location", EA="effect_allele", OA="other_allele", EAF="effect_allele_frequency", BETA="beta", SE="se", OR="odds_ratio", OR_LR="ci_lower", OR_UR="ci_upper", P="p_value", TRAIT="trait", ID="study_accession")

  }

  # apply the maps
  data.table::setnames(dat, col_map, names(col_map))
  dat[, names(dat)[!names(dat) %in% names(col_map)] := NULL]

  # lift across builds if needed
  if(!is.null(build_to) && !is.null(build_from)) {

    dat <- genepi.utils::lift(dat,
                              from=build_from,
                              to=build_to,
                              snp_col="RSID",
                              chr_col="CHR",
                              pos_col="BP",
                              ea_col="EA",
                              oa_col="OA",
                              remove_duplicates = FALSE)

  }

  # final processing steps / column additions
  dat[, nlog10P := -log10(P)]
  dat <- dat[is.numeric(BETA) &
             is.finite(BETA) &
             is.numeric(P) &
             is.finite(P) &
             is.finite(SE) &
             P < 1 &
             P > 0, ]

  # final check on minimum required standardised column names being present
  std_col_names <- c("RSID","CHR","BP","EA","OA","EAF","BETA","SE","P","N","TRAIT","ID")
  if(!all(std_col_names %in% names(dat))) {
    stop(paste0(c("Import standardisation error - minimum data columns not found. Missing:", std_col_names[which(!std_col_names %in% names(dat))]), collapse=" "))
  }

  # return standardised data
  return(dat)
}
