pathWAS_1mrscores = function(gene_PRS,
                             mr_lasso_res,
                             incl_omics = TRUE,
                             endpoint_omics = NULL
                             ) {

  PRS_iids = gene_PRS$iid
  path_qtl_ovgenes = mr_lasso_res@Exposure

  PRS_path_ovgenes = path_qtl_ovgenes[path_qtl_ovgenes %in% colnames(gene_PRS)]

  path_gene_PRS = as.data.frame(gene_PRS[, PRS_path_ovgenes])

  PRS_genedex = path_qtl_ovgenes %in% colnames(path_gene_PRS)

  for (num_expos in 1:ncol(path_gene_PRS)){

    path_gene_PRS[, num_expos] = path_gene_PRS[, num_expos] * mr_lasso_res$Estimate[PRS_genedex][num_expos]

  }

  mr_gene_PRS = data.frame(iid = PRS_iids,
                           path_gene_PRS
  )

  if (incl_omics == TRUE){

    mr_gene_PRS = merge(mr_gene_PRS, endpoint_omics, by = "iid")

  }

  return(mr_gene_PRS)

}
