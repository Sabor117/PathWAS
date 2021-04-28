pathWAS_predictR = function(predict_PRS,
                            ) {

  ovgenes_PRS = path_cohort_ovgenes[path_cohort_ovgenes %in% colnames(predict_PRS_read_rmNAs)]

  predict_PRS = as.data.frame(predict_PRS_read_rmNAs[, ovgenes_PRS])
  colnames(predict_PRS) = ovgenes_PRS

  predict_PRS_genedex = path_cohort_ovgenes %in% colnames(predict_PRS)

  for (num_expos in 1:ncol(predict_PRS)){

    predict_PRS[, num_expos] = predict_PRS[, num_expos] * mr_lasso_res$Estimate[predict_PRS_genedex][num_expos]

  }

  predict_geneno = ncol(predict_PRS)

  predict_PRS_comb = data.frame(iid = predict_PRS_read$iid, combined = rowSums(predict_PRS))

  if (length(sig_mr_genelist) > 1){

    sig_ovgenes_PRS = sig_mr_genelist[sig_mr_genelist %in% colnames(predict_PRS)]

    sig_predict_PRS = as.data.frame(predict_PRS[, sig_ovgenes_PRS])
    colnames(sig_predict_PRS) = sig_ovgenes_PRS

    sig_predict_geneno = ncol(sig_predict_PRS)

    if (ncol(sig_predict_PRS) > 0){

      colnames(sig_predict_PRS) = sig_ovgenes_PRS

      sig_predict_PRS_genedex = path_cohort_ovgenes %in% colnames(sig_predict_PRS)

      for (num_expos in 1:ncol(sig_predict_PRS)){

        sig_predict_PRS[,num_expos] = sig_predict_PRS[,num_expos] * mr_lasso_res$Estimate[sig_predict_PRS_genedex][num_expos]

      }

      sig_predict_PRS_comb = data.frame(iid = predict_PRS_read$iid, combined = rowSums(sig_predict_PRS))

    } else if (ncol(sig_predict_PRS) == 0){

      sig_predict_PRS_comb = NULL

    }

  } else if (length(sig_mr_genelist) == 0){

    sig_predict_PRS_comb = NULL
    sig_predict_geneno = 0

  }

  panel_info = olink_panPref(end_point,
                             multi_panel = FALSE)
  panel_info = data.frame(panel_info)
  panel_pref = panel_info$panel_pref[1]
  panel = panel_info$panel_name[1]

  endpoint_measures_p = fread(paste0(predictProts, "orcades_olink_proteomics_",
                                     panel,"_llod_phenotypes.tsv"),
                              data.table = FALSE)

  test_path_prot = pathway_measuR(endpoint_measures_p, end_point, panel_pref, end_point_olink$orcade_suffix)

  test_model = prs_mergeR(predict_PRS_comb, test_path_prot, end_point)
  test_sum = summary(test_model)
  test_p = lmp(test_model)

  if (!(is.null(sig_predict_PRS_comb))) {

    sig_test_model = prs_mergeR(sig_predict_PRS_comb, test_path_prot, end_point)
    sig_test_sum = summary(sig_test_model)
    sig_test_p = lmp(sig_test_model)

    outrow = data.frame(pathway = path_select, endpoint = end_point, tissue = tissue,
                        model_r2 = test_sum$r.squared, model_p = test_p, model_beta = test_sum$coefficients[2,1], gene_no = predict_geneno,
                        sigmodel_r2 = sig_test_sum$r.squared, sigmodel_p = sig_test_p, sig_beta = sig_test_sum$coefficients[2,1], sig_geneno = sig_predict_geneno)

  } else {

    outrow = data.frame(pathway = path_select, endpoint = end_point, tissue = tissue,
                        model_r2 = test_sum$r.squared, model_p = test_p, model_beta = test_sum$coefficients[2,1], gene_no = predict_geneno,
                        sigmodel_r2 = NA, sigmodel_p = NA, sig_beta = NA,  sig_geneno = sig_predict_geneno)

  }

  out_table = rbind(out_table, outrow)

}
