#' Merging clumped SNPs with the SNPs from an omics end-point
#'
#' @description
#' omics_MungeR inputs a data frame of clumped SNPs from qtl_clumpR and the summary statistics from your selected
#' end-point omics and then determines the alignment.
#'
#' @details
#' This is an intermediate step of PathWAS, required for setting up the next stages. The function inputs both data frames of the summary stats for your chosen end-point omics and the clumped SNPs,
#' most likely your output from qtl_clumpR. It combines these two data frames and produces a merged data frame
#' of the omics summary stats with an additional column which will be used to determine whether the SNP effects need to be flipped in order
#' to align the effects with those from your clumped SNPs.
#' At a minimum the input data frames must have a SNP ID column (usually rsID) and an effect allele and alternative
#' allele column. These can be defined manually as part of the function but default to "rsid,a1,a0" respectively.
#'
#' @param endpoint_omics Data frame of summary stats of your selected omics. Contains, at minimum a SNP ID column, effect allele and reference allele.
#' @param clumped_snplist Data frame of summary stats of clumped SNPs. Contains, at minimum a SNP ID column, effect allele and reference allele.
#' @param omics_SNPCol Character. Name/number of column containing the SNP IDs in the omics data frame.
#' @param omics_EffAllCol Character. Name/number of column containing the effect allele in the omics data frame.
#' @param omics_OthAllCol Character. Name/number of column containing the alternative allele in the omics data frame.
#' @param clumped_SNPCol Character. Name/number of column containing the SNP IDs in the clumped SNPs data frame.
#' @param clumped_EffAllCol Character. Name/number of column containing the effect allele in the clumped SNPs data frame.
#' @param clumped_OthAllCol Character. Name/number of column containing the alternative allele in the clumped SNPs data frame.
#'
#' @examples
#' ## Input two different data frames with different column names for merging
#' omics_MungeR(endpoint_omics = my_proteomics_stats, clumped_snplist = my_snps_dataframe, omics_SNPCol = "SNPID",
#' omics_EffAllCol = "Effect_Allele", omics_OthAllCol = "Other_Allele", clumped_SNPCol = "rsid", clumped_EffAllCol = "a1",
#' clumped_OthAllCol = "a0")
#'
omics_MungeR = function(endpoint_omics,
                        clumped_snplist,
                        omics_SNPCol = "rsid", omics_EffAllCol = "a1", omics_OthAllCol = "a0",
                        clumped_SNPCol = "rsid", clumped_EffAllCol = "a1", clumped_OthAllCol = "a0"
                        ) {

  omics_snps = endpoint_omics[endpoint_omics[, omics_SNPCol] %in% clumped_snplist,]

  omics_test = data.frame(rsid = omics_snps[, omics_SNPCol],
                          eff = omics_snps[, omics_EffAllCol],
                          ref = omics_snps[, omics_OthAllCol])

  clumped_test = data.frame(rsid = clumped_snplist[, clumped_SNPCol],
                          a1 = clumped_snplist[, clumped_EffAllCol],
                          a0 = clumped_snplist[, clumped_OthAllCol])

  ### Merge tables by rsID

  flip_table = merge(omics_test, clumped_test, by = "rsid", sort = FALSE)

  ### Utilise flip data tables function on merged table
  ### Outputs flip status as additional column

  flip_table$flip = flip.data.frame(flip_table)
  flip_table_keep = data.frame(flip_table$rsid,
                               flip_table$flip)

  colnames(flip_table_keep) = c(omics_SNPCol, "FLIP")

  omics_snps = merge(endpoint_omics,
                     flip_table_keep,
                     by = omics_SNPCol,
                     sort = FALSE)

  return(omics_snps)

}
