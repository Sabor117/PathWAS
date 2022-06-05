#' Conduct Mendelian randomisation of clumped SNPs from pathway genes on end-point omics
#'
#' @description
#'
#'
#' @details
#'
#'
#'
#' @param mr_input
#' @param end_point character. Name of end point protein (can be any string). Used for saving purposes if any save variable is set to TRUE.
#' @param path_select character. Name of pathway (can be any string). Used for saving purposes if any save variable is set to TRUE.
#'
#' @examples
#'
#' @import MendelianRandomization
#'
#' @export
pathWAS_MR = function(mr_input_check,
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
