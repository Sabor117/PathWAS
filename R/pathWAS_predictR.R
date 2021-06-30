#' Predict MR pathway models in pathway PRS
#'
#' @description
#' Creates a prediction statistic for your pathway model based on pathway PRS
#'
#' @details
#' pathWAS_predictR is the final step in PathWAS. It requires the output from pathWAS_MR which contains the MR-derived
#' gene-exposure values. As well as this it requires two important aspects for creating the final PathWAS model.
#' 1. It requires your end-point proteomics measures. The end-point selected at the start of PathWAS must be input
#' along with a data frame of iids + omics measures. This MUST be in the format of 2 columns: iid, metabolite/protein (where the name
#' of the second column is in the format of: "XXX_omic" where "XXX" can be anything).
#' 2. It requires polygenic risk scores (PRS) for every gene in your pathway. This has to be created in the same
#' cohort for which you have your omics measurements. I.e. you must take the QTLs used for the previous
#' stages of PathWAS and then use those SNPs to create weights using the genotype of your cohort.
#' For this segment we recommend the snp_ldpred2_auto() function from LDpred2 to create polygenic weights:
#' https://privefl.github.io/bigsnpr/articles/LDpred2.html. Then we recommend using PRsice2 to create PRS from
#' each polygenic weight.
#' In theory you should then have a data frame of PRS for each gene in your pathway in the format of:
#' iid, GENE1_PRS, GENE2_PRS, GENE3_PRS (where the 1st column MUST be "iid" and every subsequent column must be
#' the name of a gene).
#' You will also need to input the list of genes for your pathway, in the same format as the column names from your
#' data frame of PRS.
#' Lastly, there is the option of running this on all of the genes from your pathway or on those genes which
#' were significant within the MR. By setting the run_sig_MR argument to TRUE you will output a model for all genes
#' within the pathway that are available and also for those which were significant in the MR analysis, providing
#' you with two seprate models. If this is set to true, you will also need to input the sig_mr_genelist which
#' is output from the pathWAS_MR.
#'
#' @param predict_PRS data frame. 1st column must be titled "iid". Every subsequent column is a PRS for a single gene.
#' @param path_qtl_ovgenes list. List of all genes overlapping between your pathway and the QTLs available. These names must be in the same format as the column names for your PRS.
#' @param path_select character. Name of the pathway for the analysis.
#' @param mr_lasso_res LASSO model. Output [1] of the pathWAS_MR function.
#' @param endpoint_omics data frame. 2 columns. The first MUST be titled "iid" and the second is the omics measurement of your end-point. These IIDs must be the same as for the PRS data frame.
#' @param end_point character. Name of the omics end-point measured and used as the proxy for pathway functionality.
#' @param run_sig_MR logical. Also run a prediction on only significant exposures (genes) from the MR. Default is FALSE.
#' @param sig_mr_genelist list. If run_sig_MR == TRUE, you must include this. Output [2] from the pathWAS_MR function.
#'
#' @export
pathWAS_predictR = function(predict_PRS,
                            path_qtl_ovgenes,
                            path_select,
                            mr_lasso_res,
                            endpoint_omics,
                            end_point,
                            run_sig_MR = FALSE,
                            sig_mr_genelist = NULL
                            ) {

  PRS_iids = predict_PRS$iid

  PRS_path_ovgenes = path_qtl_ovgenes[path_qtl_ovgenes %in% colnames(predict_PRS)]

  all_predict_PRS = as.data.frame(predict_PRS[, PRS_path_ovgenes])

  predict_PRS_genedex = path_qtl_ovgenes %in% colnames(all_predict_PRS)

  for (num_expos in 1:ncol(all_predict_PRS)){

    all_predict_PRS[, num_expos] = all_predict_PRS[, num_expos] * mr_lasso_res@Estimate[predict_PRS_genedex][num_expos]

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

          sig_predict_PRS[,num_expos] = sig_predict_PRS[,num_expos] * mr_lasso_res@Estimate[sig_predict_PRS_genedex][num_expos]

        }

        sig_predict_PRS_comb = data.frame(iid = PRS_iids,
                                          combined = rowSums(sig_predict_PRS))

      } else if (ncol(sig_predict_PRS) == 0){

        sig_predict_PRS_comb = NULL

      }

    } else if (length(sig_mr_genelist) <= 1){

      sig_predict_PRS_comb = NULL
      sig_predict_geneno = length(sig_mr_genelist)

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
