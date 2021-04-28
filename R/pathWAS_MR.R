pathWAS_MR = function(path_snplist,
                         clumped_snps,
                         qtl_sumstats,
                         geneCol = NULL,
                         omics_snps,
                         omics_SNPCol= "rsid", omics_BetaCol = "beta1", omics_SECol = "se",
                         end_point = NULL, path_select = NULL,
                         save_MRInput = FALSE,
                         save_MRInLoc = NULL,
                         save_MROutput = FALSE,
                         save_MROutLoc = NULL,
                         save_MRExps = FALSE,
                         save_MRExpsLoc = NULL
                        ) {

  clumped_snplist = clumped_snps$rsid
  path_cohort_ovgenes = unique(path_snplist$gene)

  snp_beta_matrix = data.frame(matrix(vector(), 0, length(path_cohort_ovgenes),
                                      dimnames = list(c(), path_cohort_ovgenes)),
                               stringsAsFactors = F)

  snp_se_matrix = data.frame(matrix(vector(), 0, length(path_cohort_ovgenes),
                                    dimnames = list(c(), path_cohort_ovgenes)),
                             stringsAsFactors = F)

  for (nsnp in 1:length(clumped_snplist)){

    currSnp = clumped_snplist[nsnp]

    snpRow_beta = c()
    snpRow_se = c()

    for (ngene in 1:length(path_cohort_ovgenes)){

      currGene = path_cohort_ovgenes[ngene]

      ### EDIT: Individual gene file

      if (grepl("%%%", qtl_sumstats)){

        all_gene_sumstats = fread(gsub("%%%", currGene, qtl_sumstats),
                                  data.table = FALSE)

      } else {

        all_qtl_sumstats = fread(qtl_sumstats,
                                  data.table = FALSE)

        all_gene_sumstats = all_qtl_sumstats[currGene %in% all_qtl_sumstats[, geneCol],]

      }

      snp_gene_beta = all_gene_sumstats$beta1[all_gene_sumstats$rsid == currSnp]
      snp_gene_se = all_gene_sumstats$beta1[all_gene_sumstats$rsid == currSnp]

      if (length(snp_gene_beta) == 0){

        snp_gene_beta = 0.0000001
        snp_gene_se = 1

      }

      snpRow_beta = c(snpRow_beta, snp_gene_beta)
      snpRow_se = c(snpRow_se, snp_gene_se)

    }

    snpRow_beta = t(data.frame(snpRow_beta))
    colnames(snpRow_beta) = path_cohort_ovgenes
    rownames(snpRow_beta) = currSnp

    snp_beta_matrix = rbind(snp_beta_matrix, snpRow_beta)

    snpRow_se = t(data.frame(snpRow_se))
    colnames(snpRow_se) = path_cohort_ovgenes
    rownames(snpRow_se) = currSnp

    snp_se_matrix = rbind(snp_se_matrix, snpRow_se)

  }

  snp_beta_matrix = as.matrix(snp_beta_matrix)
  snp_se_matrix = as.matrix(snp_se_matrix)

  if (ncol(snp_beta_matrix) == 1){

    matrix_cols = colnames(snp_beta_matrix)

    snp_beta_matrix = as.matrix(data.frame(snp_beta_matrix[omics_snps[, omics_SNPCol],]))
    snp_se_matrix = as.matrix(data.frame(snp_se_matrix[omics_snps[, omics_SNPCol],]))

    colnames(snp_beta_matrix) = matrix_cols
    colnames(snp_se_matrix) = matrix_cols

  } else {

    omics_snps_list = omics_snps[, omics_SNPCol]

    snp_beta_matrix = snp_beta_matrix[omics_snps[, omics_SNPCol],]
    snp_se_matrix = snp_se_matrix[omics_snps[, omics_SNPCol],]

  }

  omics_betas = omics_snps[, omics_BetaCol] * omics_snps$FLIP
  omics_se = omics_snps[, omics_SECol]

  heading("Matrices made. Creating MR input.")

  mr_input = mr_mvinput(bx = snp_beta_matrix,
                        bxse = snp_se_matrix,
                        by = omics_betas,
                        byse = omics_se)

  heading("MR input made. Running MR.")

  if (save_MRInput == TRUE) {

    saveRDS(mr_input,
            file = paste0(save_MRInLoc, end_point, "_", path_select, "_MR_input.rds"))


  }

  mr_lasso_res = mr_mvlasso(mr_input)

  mr_lasso_res@Exposure = colnames(snp_beta_matrix)

  if (save_MROutput == TRUE) {

    saveRDS(mr_lasso_res,
          file = paste0(save_MROutLoc, end_point, "_", path_select, "_MR_output.rds"))

  }

  heading("MR complete.")

  mr_lasso_names = list(genes = colnames(snp_beta_matrix), rsids = rownames(snp_beta_matrix))

  if (save_MRExps == TRUE) {

    saveRDS(mr_lasso_names,
          file = paste0(save_MRExpsLoc, end_point, "_", path_select, "_MR_exps.rds"))

  }

  mr_sig_index = mr_lasso_res$Pvalue < 0.05
  sig_mr_genelist = path_cohort_ovgenes[mr_sig_index]

  return(c(mr_lasso_res, sig_mr_genelist))

}
