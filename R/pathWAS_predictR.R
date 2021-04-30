pathWAS_predictR = function(predict_PRS,
                            path_qtl_ovgenes, path_select,
                            mr_lasso_res,
                            endpoint_omics, end_point,
                            run_sig_MR = FALSE,
                            sig_mr_genelist = NULL
                            ) {

  PRS_iids = predict_PRS$iid

  PRS_path_ovgenes = path_qtl_ovgenes[path_qtl_ovgenes %in% colnames(predict_PRS)]

  all_predict_PRS = as.data.frame(predict_PRS[, PRS_path_ovgenes])
#  colnames(predict_PRS) = PRS_path_ovgenes

  predict_PRS_genedex = path_qtl_ovgenes %in% colnames(all_predict_PRS)

  for (num_expos in 1:ncol(all_predict_PRS)){

    all_predict_PRS[, num_expos] = all_predict_PRS[, num_expos] * mr_lasso_res$Estimate[predict_PRS_genedex][num_expos]

  }

  predict_geneno = ncol(all_predict_PRS)

  predict_PRS_comb = data.frame(iid = PRS_iids,
                                combined = rowSums(all_predict_PRS))

  if (run_sig_MR == TRUE){

    if (length(sig_mr_genelist) > 1){

      sig_ovgenes_PRS = sig_mr_genelist[sig_mr_genelist %in% colnames(predict_PRS)]

      sig_predict_PRS = as.data.frame(predict_PRS[, sig_ovgenes_PRS])

      sig_predict_geneno = ncol(sig_predict_PRS)

      if (ncol(sig_predict_PRS) > 0){

        colnames(sig_predict_PRS) = sig_ovgenes_PRS

        sig_predict_PRS_genedex = path_cohort_ovgenes %in% colnames(sig_predict_PRS)

        for (num_expos in 1:ncol(sig_predict_PRS)){

          sig_predict_PRS[,num_expos] = sig_predict_PRS[,num_expos] * mr_lasso_res$Estimate[sig_predict_PRS_genedex][num_expos]

        }

        sig_predict_PRS_comb = data.frame(iid = PRS_iids,
                                          combined = rowSums(sig_predict_PRS))

      } else if (ncol(sig_predict_PRS) == 0){

        sig_predict_PRS_comb = NULL

      }

    } else if (length(sig_mr_genelist) == 0){

      sig_predict_PRS_comb = NULL
      sig_predict_geneno = 0

    }

  }

  test_model = prs_mergeR(predict_PRS_comb,
                          endpoint_omics,
                          end_point)
  test_sum = summary(test_model)
  test_p = lmp(test_model)

  if (run_sig_MR == FALSE){

    outrow = data.frame(pathway = path_select, endpoint = end_point,
                        model_r2 = test_sum$r.squared, model_p = test_p, model_beta = test_sum$coefficients[2,1], gene_no = predict_geneno)

  } else if (!(is.null(sig_predict_PRS_comb))) {

    sig_test_model = prs_mergeR(sig_predict_PRS_comb,
                                endpoint_omics,
                                end_point)
    sig_test_sum = summary(sig_test_model)
    sig_test_p = lmp(sig_test_model)

    outrow = data.frame(pathway = path_select, endpoint = end_point,
                        model_r2 = test_sum$r.squared, model_p = test_p, model_beta = test_sum$coefficients[2,1], gene_no = predict_geneno,
                        sigmodel_r2 = sig_test_sum$r.squared, sigmodel_p = sig_test_p, sig_beta = sig_test_sum$coefficients[2,1], sig_geneno = sig_predict_geneno)

  } else {

    outrow = data.frame(pathway = path_select, endpoint = end_point,
                        model_r2 = test_sum$r.squared, model_p = test_p, model_beta = test_sum$coefficients[2,1], gene_no = predict_geneno,
                        sigmodel_r2 = NA, sigmodel_p = NA, sig_beta = NA,  sig_geneno = sig_predict_geneno)

  }

  return(outrow)

}
