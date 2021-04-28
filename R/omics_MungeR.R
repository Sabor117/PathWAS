omics_MungeR = function(endpoint_omics,
                        clumped_snplist,
                        omics_SNPCol = "rsid", omics_EffAllCol = "a1", omics_OthAllCol = "a0") {

  omics_snps = endpoint_omics[endpoint_omics[, omics_SNPCol] %in% clumped_snplist,]

  omics_test = data.frame(rsid = omics_snps[, omics_SNPCol],
                          eff = omics_snps[, omics_EffAllCol],
                          ref = omics_snps[, omics_OthAllCol])

  clumped_test = data.frame(rsid = clumped_snps$rsid,
                          a1 = clumped_snps$a1,
                          a0 = clumped_snps$a0)

  ### Merge tables by rsID

  flip_table = merge(omics_test, clumped_test, by = "rsid", sort = FALSE)

  ### Utilise flip data tables function on merged table
  ### Outputs flip status as additional column

  flip_table$flip = flip.data.frame(flip_table)
  flip_table_keep = data.frame(flip_table$rsid,
                               flip_table$flip)

  colnames(flip_table_keep) = c(omics_SNPCol, "FLIP")

  omics_snps = merge(omics_snps,
                     flip_table_keep,
                     by = omics_SNPCol,
                     sort = FALSE)

  return(omics_snps)

}
