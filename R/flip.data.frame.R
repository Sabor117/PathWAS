#' Check and flip the alleles of a data frame
#'
#' @description
#' Lower level function which takes an input data-frame with alleles from two separate genotypes
#' and then compares the alleles and flips them if necessary.
#'
#' @details
#' flip.data.frame is called within other major functions and is required for aligning the direction of
#' effects between two different datasets for performing a regression between them. The function (and the
#' other functions it calls) operates in two ways: first checking whether alleles a1/a0 in one dataset match
#' a1/a0 in the other. And then it switches the strand of the alleles being checked to make sure that these do
#' not mismatch with those on the other side either.
#'
#' @param df data frame. With 5 columns in order: rsID, a1, a0 (genotype to be checked), a1, a0 (genotype to be checked against)
#'
#' @export
flip.data.frame = function(df) {

  apply(df[2:5], 1, function(x) flip.allele(x[1], x[2], x[3], x[4]))

}
