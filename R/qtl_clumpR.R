#' Clump QTL SNPs based on gene
#'
#' @description qtl_clumpR returns two outputs: a compiled data frame of SNPs for a given pathway along with selected clumped SNPs
#'
#' @details This function takes as an input several large data frames. You must input a data frame of your QTLs
#' (e.g. qQTLs/pQTLs) and also a list of your genes for your pathway.
#' This can optionally be filtered to only include significant QTLs/eGenes (and we recommemnd this).
#'
#' The QTLs must (for this function) at minimum include 4 columns:
#' rsid column - e.g. rs1234.
#' snpid column - a unique identifier for every SNP. We suggest the format "chromosome:position_a1_a2" but this is
#' not compulsory.
#' p - the QTL p-value for each SNP.
#' gene - each SNP must be labeled according to which gene it is a QTL for. Many eQTL data sets use Ensembl IDs as
#' a default. The name of this column can be chosen as a variable of the function.
#'
#' SNPs should be included for every gene available in the pathway. Where SNPs overlap between genes, this
#' function will select said SNP only once based on the lowest QTL p-value.
#'
#' You must also include a BioMart file (from https://www.ensembl.org/biomart/martview/ or using package biomaRt)
#' in order to convert the gene names between the names used in the QTLs and Entrez (KEGG) gene IDs.
#'
#' The clumping is performed using the package ieugwasr (https://rdrr.io/github/MRCIEU/ieugwasr/).
#' This requires access to a locally installed version of plink and access to a reference genotype
#' in plink format (I.e. bed/bim/fam).
#'
#' For the requirement of a reference genotype you should input the location of the files + the prefix of your
#' plink files (e.g. if you have /opt/data/reference/cohort_new_allchr.bed,
#' /opt/data/reference/cohort_new_allchr.bim, cohort_new_allchr.fam use
#' bfile = "/opt/data/reference/cohort_new_allchr"). If your plink files are divided by chromosome
#' (e.g. cohort_new_chr1.bed : cohort_new_chr2.bed) then input the prefix with "%%%" replacing the
#' chromosome number: bfile = "/opt/data/reference/cohort_new_chr%%%"
#'
#' This function will output (as a list) two things: the filtered data frame of all the QTL SNPs for the pathway
#' and a data frame of the clumped SNPs (usually only a few per locus).
#'
#' @param end_point character. Name of gene in HGNC format (or equivalent). Required for excluding the gene from
#' the list of QTLs used.
#' @param path_select character. Name of pathway. Format is unimportant, only used for aesthetics purposes.
#' @param path_gene_list list. The output of either smple_paths or genepath_ListR. A list of Entrez IDs of all genes in your pathway.
#' @param biomart_map data frame. Data frame of gene names including at least entrezgene_id, external_gene_name
#' and whatever format of gene name is used in the QTLs (e.g. Ensembl ID).
#' @param all_snps data frame. All SNPs for QTLs to be filtered in analysis. May be significant QTLs only.
#' @param all_snps_genecol character or numeric. Name of the column containing the gene name from the QTLs in all_snps. Default is
#' "gene_ensembl" as many QTL datasets tend to use Ensembl IDs of some description. Can also be the number of the
#' column containing the gene names.
#' @param qtl_gene_identifier character or numeric. Equivalent naming format from the SNPs in the BioMart map. E.g.
#' if the QTLs use UniProt or Ensembl IDs, then this should be the equivalent from BioMart.
#' @param bfile character. Location of the Plink files of your reference genotype.
#' @param plink_bin character. Location of the local version of the plink executable.
#'
#' @import data.table ieugwasr
#'
#' @export
qtl_clumpR = function(end_point, path_select,
                      path_gene_list,
                      biomart_map,
                      all_snps,
                      all_snps_genecol = "gene_ensembl",
                      qtl_gene_identifier = "ensembl_gene_id",
                      bfile,
                      plink_bin,
                      MAF_filter = NA,
                      MAF_col = "MAF",
                      def_tmpDir = tempdir()
) {

  biomart_map = data.table::fread(biomart_map,
                                  data.table = FALSE)

  genelist_frame = data.frame(unique(biomart_map[biomart_map$entrezgene_id %in% path_gene_list,]))
  genelist_frame = genelist_frame[!(genelist_frame$external_gene_name == end_point),]

  path_genes_name = unique(genelist_frame$external_gene_name)
  path_genes_identifier = unique(genelist_frame[,qtl_gene_identifier])

  path_snplist = all_snps[all_snps[,all_snps_genecol] %in% path_genes_identifier,]

  if (!(is.na(MAF_filter))){

    if (!(MAF_col %in% colnames(path_snplist))){

      stop("Please enter MAF column name if you wish to filter by MAF.")

    }

    path_snplist = path_snplist[path_snplist[,MAF_col] > MAF_filter,]

  }

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

    clumped_snps = pathwas_ld_clump(dat = path_snplist,
                                      bfile = bfile,
                                      plink_bin = plink_bin,
                                      def_tmpDir = def_tmpDir
    )

  } else {

    for (numchrom in 1:length(unique(path_snplist$chr))){

      chrom = unique(path_snplist$chr)[numchrom]

      path_snplist_chrom = data.frame(path_snplist[path_snplist$chr == chrom,])
      colnames(path_snplist_chrom)[colnames(path_snplist_chrom) == "p"] = "pval"

      curr_clumped_snps = pathwas_ld_clump(dat = path_snplist_chrom,
                                             bfile = gsub("%%%", chrom, bfile),
                                             plink_bin = plink_bin,
                                             def_tmpDir = def_tmpDir
      )

      clumped_snps = rbind(clumped_snps, curr_clumped_snps)

    }
  }

  return(clumped_snps)
}
