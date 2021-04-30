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
