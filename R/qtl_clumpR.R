qtl_clumpR = function(end_point, path_select,
                        path_gene_list,
                        biomart_map,
                        all_snps,
                        all_snps_genecol = "gene_ensembl",
                        qtl_gene_identifier = "ensembl_gene_id",
                        tissue = "Complete",
                        bfile,
                        plink_bin
                        ) {

  require(ieugwasr)
  require(stringr)

  tiss_genes_kegg = path_gene_list$gene[path_gene_list$tissue == tissue]
  tiss_genes_entrez = str_split_fixed(tiss_genes_kegg, ":", 2)[,2]

  genelist_frame = data.frame(unique(biomart_map[biomart_map$entrezgene_id %in% tiss_genes_entrez,]))
  genelist_frame = genelist_frame[!(genelist_frame$external_gene_name == end_point),]

  tiss_genes_name = unique(genelist_frame$external_gene_name)
  tiss_genes_identifier = unique(genelist_frame[,qtl_gene_identifier])

  path_snplist = all_snps[all_snps[,all_snps_genecol] %in% tiss_genes_identifier,]

  repeat_snps = path_snplist$snpid[duplicated(path_snplist$snpid)]
  repeat_snps_info = path_snplist[path_snplist$snpid %in% repeat_snps,]
  rem_snps = path_snplist[!(path_snplist$snpid %in% repeat_snps),]

  lowest_p_snps = data.frame(setDT(repeat_snps_info)[, .SD[which.min(p)], by = snpid])

  path_snplist = rbind(rem_snps, lowest_p_snps)

  if (nrow(path_snplist) == 0){

    heading(paste0("No SNPs present for pathway: ", path_select, " with end point: ", end_point))
    stop("Nae SNPs.")

  }

  clumped_snps_cols = colnames(path_snplist)
  clumped_snps_cols[clumped_snps_cols == "p"] = "pval"
  clumped_snps_cols[length(clumped_snps_cols) + 1] = "id"

  clumped_snps = data.frame(matrix(vector(), 0, length(clumped_snps_cols),
                            dimnames = list(c(), clumped_snps_cols)),
                            stringsAsFactors = F)

  if (!(grepl("%%%", bfile))){

    clumped_snps = ieugwasr::ld_clump(dat = path_snplist,
                                           bfile = bfile,
                                           plink_bin = plink_bin
                                            )

  } else {

    for (numchrom in 1:length(unique(path_snplist$chr))){

      chrom = unique(path_snplist$chr)[numchrom]

      path_snplist_chrom = data.frame(path_snplist[path_snplist$chr == chrom,])
      colnames(path_snplist_chrom)[colnames(path_snplist_chrom) == "p"] = "pval"

      curr_clumped_snps = ieugwasr::ld_clump(dat = path_snplist_chrom,
                                             bfile = gsub("%%%", chrom, bfile),
                                             plink_bin = plink_bin
      )

      clumped_snps = rbind(clumped_snps, curr_clumped_snps)

    }
  }

  return(c(path_snplist, clumped_snps))
}
