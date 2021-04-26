qtl_clumpR = function(end_point, path_select,
                        path_gene_list,
                        biomart_map,
                        all_snps,
                        snps_gene_col,
                        qtl_gene_identifier = "ensembl_gene_id",
                        tissue = "Complete"
                          ) {

  require(ieugwasr)

  tiss_genes_kegg = path_gene_list$gene[path_gene_list$tissue == tissue]
  tiss_genes_entrez = str_split_fixed(tiss_genes_kegg, ":", 2)[,2]

  genelist_frame = data.frame(unique(biomart_map[biomart_map$entrezgene_id %in% tiss_genes_entrez,]))
  genelist_frame = genelist_frame[!(genelist_frame$external_gene_name == end_point),]

  tiss_genes_name = unique(genelist_frame$external_gene_name)
  tiss_genes_identifier = unique(genelist_frame[qtl_gene_identifier,])

  path_snplist = all_snps[all_snps[qtl_gene_identifier,] %in% tiss_genes_identifier,]

  repeat_snps = path_snplist$snpid[duplicated(path_snplist$snpid)]
  repeat_snps_info = path_snplist[path_snplist$snpid %in% repeat_snps,]
  signif_rem = path_snplist[!(path_snplist$snpid %in% repeat_snps),]

  lowest_p_snps = data.frame(setDT(repeat_snps_info)[, .SD[which.min(p)], by = snpid])

  path_snplist = rbind(signif_rem, lowest_p_snps)

  if (nrow(path_sig_snplist) == 0){

    heading(paste0("No SNPs present for pathway: ", path_select, " with end point: ", end_point))
    stop("Nae SNPs.")

  }

  clumped_snps_cols = colnames(path_snplist)
  clumped_snps_cols[11] = "pval"
  clumped_snps_cols[14] = "id"

  clumped_snps = data.frame(matrix(vector(), 0, 14,
                            dimnames = list(c(), clumped_snps_cols)),
                            stringsAsFactors = F)

  for (numchrom in 1:length(unique(path_snplist$chr))){

    chrom = unique(path_sig_snplist$chr)[numchrom]

    path_sig_snplist_reform = data.frame(path_sig_snplist[path_sig_snplist$chr == chrom,])
    colnames(path_sig_snplist_reform)[11] = "pval"

    curr_clumped_snps = ieugwasr::ld_clump(dat = path_sig_snplist_reform,
                                           bfile = paste0(refPlink, "ukbb_chr", chrom, "_10000_random_unrelated_white_british"),
                                           plink_bin = "/exports/igmm/software/pkg/el7/apps/plink/1.90b4/plink"
                                            )

    clumped_snps = rbind(clumped_snps, curr_clumped_snps)

  }

  clumped_snplist = clumped_snps$rsid
  path_cohort_ovgenes = unique(path_sig_snplist$gene)

  if (length(clumped_snplist) <= 1){

    heading(paste0("Only one SNP present after clumping for pathway: ", path_select, " with end point: ", end_point))
    stop("Nae enough SNPs.")

  }
}
