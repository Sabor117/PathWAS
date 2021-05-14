#' Creates model for pathway regressed against end-point omics
#'
#' @description
#' prs_mergeR is a lower level function for pathWAS_predictR which merges a data frame of PRS for a given pathway
#' and a data frame of selected omics (and the name of the end-point) and then performs a regression on this merged data frame.
#'
#' @details
#' In pathWAS_predictR this function is used to merge two data frames. One data frame is 2 columns: IID and pathway PRS
#' (where pathway PRS is a summed PRS for every gene PRS in the pathway). The second data frame is also 2 columns:
#' IID and the proteomics measurements for your selected end-point.
#' The function will merge these two data frames by IID and then regresses the omics values by the summed PRS.
#' To function correctly the IIDs must overlap between both data frames (I.e. PRS must be made in the same individuals
#' for which you have your end-point omics).
#'
#' @param pathway_scores data frame. At LEAST 2 columns: MUST contain "iid" column. Every other column is either individual gene PRS or a summed pathway PRS in one column
#' @param omics data frame. ONLY 2 columns. MUST contain "iid" column overlapping with pathway_scores. Other column is the end-point omics.
#' @param gene character. Name of the end-point protein used in the analysis.
#'
prs_mergeR = function(pathway_scores, omics, gene){

  ### Combine scores and proteomics

  test_path_frame = merge(omics, pathway_scores, by = "iid")
  test_path_frame = test_path_frame[, apply(test_path_frame, 2, function(x) !any(is.na(x)))] ### Remove NAs

  ### Define model forumlae based on protein name
  ### Then create model and produce summary statistics

  model_formulae = as.formula(paste0(gene, "_omic ~ .")) ### Output relies on input from previous package functions

  model = lm(model_formulae, data = test_path_frame[,-1])

  return(model)

}
