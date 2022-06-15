#' Removes invalid SNPs from MR input
#'
#' @description
#' Takes in MR input created by pathWAS_MR and checks which (if any) SNPs are deemed invalid by the MendelianRandomization
#' package,
#'
#' @details
#' The MendelianRandomization::mr_mvlasso() function used in the PathWAS package produces an output with a table
#' of exposures (genes) and estimates (weights/betas) based on SNPs input into the function.
#' Occasionally, this output will have a different length of estimates and exposures due to removal of SNPs deemed
#' as "invalid" by the function. This causes issues with downstream analyses as it means that some weights
#' are missing from the output table, and thus cannot be used in creation of PRS. It also means that reading/extracting
#' information from the table can cause errors in code. As the rationable behind what SNPs are deemed "invalid" is slightly arcane,
#' this function is designed to use the MR input object used in mr_mvlasso() and find the SNPs which cause
#' the erroneous output, so that the mr_mvlasso() function can be run without those SNPs in the first place.
#'
#' @param mr_input MRInput object. Object created by the MendelianRandomization::mr_mvinput() function
#' @param end_point character. Name of end point protein (can be any string). Used for saving purposes if any save variable is set to TRUE.
#' @param path_select character. Name of pathway (can be any string). Used for saving purposes if any save variable is set to TRUE.
#'
#' @examples
#'
#' @import MendelianRandomization
#'
#' @export
pathWAS_snpCheck = function(mr_input_check,
                        verbose = TRUE
                        ) {

  if (verbose == TRUE){

    heading("We found a bad SNP! May we burn her!")

  }

  snplist = rownames(mr_input_check@betaX)

  mr_lasso_res = MendelianRandomization::mr_mvlasso(mr_input_check)

  invalidSNPs = rownames(mr_input_check@betaX)[!(mr_input_check$snps %in% mr_lasso_res@ValidSNPs)]

  return(invalidSNPs)

}
