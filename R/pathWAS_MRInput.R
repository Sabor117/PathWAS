#' Create MR Input object for MVMR
#'
#' @description
#' Inputs and combines aspects from your QTL SNPs, clumped SNPs, pathway genes, and omics summary
#' stats (plus flip status) and creates MR input object.
#'
#' @details
#' This is a shortened version of the pathWAS_MR() function designed to solely output the same MR input
#' object used in the mr_mvinput() function if PathWAS. Created primarily for use with the
#' pathWAS_snpCheck function in mind.
#'
#' @param genelist list. List of genes extracted from database for your pathway which overlap with your QTL data.
#' @param clumped_snps data frame. The data frame of clumped SNPs output from qtl_clumpR.
#' @param qtl_sumstats character. File location for the QTLs used (either divided by gene or complete data).
#' @param geneCol character. Name or number of column containing the name of the gene for each SNP in the QTL data. Default is "gene". Ensure the format of the name is the same as your gene list.
#' @param omics_snps data frame. Summary stats of end-point omcis, now including FLIP column frm omics_MungeR
#' @param omics_SNPCol character. Name/number of column containing the SNP ID in the omics data frame.
#' @param omics_BetaCol character. Name/number of column containing the effect size in the omics data frame.
#' @param omics_SECol character. Name/number of column containing the standard error in the omics data frame.
#'
#' @examples
#' ## Having created a list of genes, clumped SNPs and munged your Omics SNPs you can input them in the following way, along with adding the option to save the MR output:
#' pathWAS_MR(genelist, clumped_snps, qtl_sumstats = "/opt/localdir/gene_qtls/file1_sumstats_$$$.tsv", geneCol = "gene_ensembl", omics_snps, save_MROutput = TRUE, save_MROutLoc = "/opt/localdir/mr_outputs/", path_select = "nod_signalling", end_point = "IL18")
#'
#' @import data.table MendelianRandomization
#'
#' @export
#'

pathWAS_MRInput = function(genelist,
                      clumped_snps,
                      qtl_sumstats,
                      geneCol = "gene",
                      omics_snps,
                      omics_SNPCol= "rsid", omics_BetaCol = "beta1", omics_SECol = "se",
                      end_point = "endprotein",
                      verbose = TRUE
                      ) {

  clumped_snplist = unique(clumped_snps$rsid)
  path_cohort_ovgenes = unique(genelist[!(genelist %in% end_point)])

  snp_beta_matrix = data.frame(matrix(ncol = (length(path_cohort_ovgenes)), nrow = 0))
  colnames(snp_beta_matrix) = c(path_cohort_ovgenes)

  snp_se_matrix = data.frame(matrix(ncol = (length(path_cohort_ovgenes)), nrow = 0))
  colnames(snp_se_matrix) = c(path_cohort_ovgenes)

  if (verbose == TRUE){

    heading("Creating SNP beta and SE matrices.")

  }

  for (nsnp in 1:length(clumped_snplist)){

    if (nsnp %% 100 == 0 && verbose == TRUE){

      cat(paste0("\nNow working on SNP: ", nsnp, "\n====\n"))

    }

    currSnp = clumped_snplist[nsnp]

    snpRow_beta = c()
    snpRow_se = c()

    for (ngene in 1:length(path_cohort_ovgenes)){

      all_gene_sumstats = NULL

      currGene = path_cohort_ovgenes[ngene]

      ### EDIT: Individual gene file

      if (grepl("%%%", qtl_sumstats)){

        if (file.exists(gsub("%%%", currGene, qtl_sumstats))){

          all_gene_sumstats = data.table::fread(gsub("%%%", currGene, qtl_sumstats),
                                                data.table = FALSE)

        } else {

          snp_gene_beta = NULL

        }
      } else {

        all_qtl_sumstats = data.table::fread(qtl_sumstats,
                                             data.table = FALSE)

        all_gene_sumstats = all_qtl_sumstats[currGene %in% all_qtl_sumstats[, geneCol],]

        if (nrow(all_gene_sumstats) == 0){

          snp_gene_beta = NULL

        }
      }

      if (exists("all_gene_sumstats")){

        snp_gene_beta = all_gene_sumstats$beta1[all_gene_sumstats$rsid == currSnp]
        snp_gene_se = all_gene_sumstats$se[all_gene_sumstats$rsid == currSnp]

      }

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

  if (all(snp_beta_matrix == 0.0000001)){

    stop("No SNPs overlap between QTLs and omics.\n\n========\n")

  }

  ### Checking for empty rows

  row_rms = c()

  for (rowcheck in 1:nrow(snp_se_matrix)){

    if (all(as.vector(snp_se_matrix[rowcheck, ]) == 1)){

      row_rms = c(row_rms, as.numeric(rowcheck))

      if (verbose == TRUE){

        cat(paste0("\n\nSNP ", rownames(snp_se_matrix)[rowcheck], " not found related to any genes.\n___________\n\n"))

      }
    }
  }

  if (length(row_rms) > 0){

    keep_rows = rownames(snp_beta_matrix)[-row_rms]

    snp_beta_matrix = as.matrix(snp_beta_matrix[-row_rms, ])
    snp_se_matrix = as.matrix(snp_se_matrix[-row_rms, ])

    rownames(snp_beta_matrix) = keep_rows
    rownames(snp_se_matrix) = keep_rows

  }

  ### Checking for empty columns

  col_rms = c()

  for (colcheck in 1:ncol(snp_se_matrix)){

    #cat(paste0(matrixStats::count(!(as.vector(snp_se_matrix[, colcheck]) == 1), TRUE), " SNP values.\n"))

    if (all(as.vector(snp_se_matrix[, colcheck]) == 1)){

      col_rms = c(col_rms, as.numeric(colcheck))

      if (verbose == TRUE){

        cat(paste0("\n\nGene ", colnames(snp_se_matrix)[colcheck], " has no overlapping SNPs.\n___________\n\n"))

      }
    }
  }

  if (length(col_rms) > 0){

    keep_cols = colnames(snp_beta_matrix)[-col_rms]

    snp_beta_matrix = as.matrix(snp_beta_matrix[, -col_rms])
    snp_se_matrix = as.matrix(snp_se_matrix[, -col_rms])

    colnames(snp_beta_matrix) = keep_cols
    colnames(snp_se_matrix) = keep_cols

  }

  if (ncol(snp_beta_matrix) == 1){

    warning("Only 1 gene remaining in pathway. Unlikely to be good model.")

    matrix_col = colnames(snp_beta_matrix)

    snp_beta_matrix = as.matrix(data.frame(snp_beta_matrix[rownames(snp_beta_matrix) %in% omics_snps[, omics_SNPCol],]))
    snp_se_matrix = as.matrix(data.frame(snp_se_matrix[rownames(snp_se_matrix) %in% omics_snps[, omics_SNPCol],]))

    colnames(snp_beta_matrix) = matrix_col
    colnames(snp_se_matrix) = matrix_col

    if (nrow(snp_beta_matrix) < 3){

      cat(paste0("QTL SNP to beta matrix:\n\n"))
      print(snp_beta_matrix)
      cat("\n\n")
      stop("Too few SNPs to use for MR.\n=====\n\n")

    }

  } else {

    omics_snps_list = omics_snps[, omics_SNPCol]

    snp_beta_matrix = snp_beta_matrix[rownames(snp_beta_matrix) %in% omics_snps[, omics_SNPCol],]
    snp_se_matrix = snp_se_matrix[rownames(snp_se_matrix) %in% omics_snps[, omics_SNPCol],]

  }

  ### Putting SNPs in same order

  snp_beta_matrix = snp_beta_matrix[order(rownames(snp_beta_matrix)),]
  snp_se_matrix = snp_se_matrix[order(rownames(snp_se_matrix)),]

  omics_snps = omics_snps[omics_snps[, omics_SNPCol] %in% rownames(snp_beta_matrix),]
  omics_snps = omics_snps[order(omics_snps[, omics_SNPCol]),]

  if (verbose == TRUE){

    heading("Matrices made. Creating MR input.")

  }

  omics_betas = omics_snps[, omics_BetaCol] * omics_snps$FLIP
  omics_se = omics_snps[, omics_SECol]

  mr_input = MendelianRandomization::mr_mvinput(bx = snp_beta_matrix,
                                                bxse = snp_se_matrix,
                                                by = omics_betas,
                                                byse = omics_se)

  return(mr_input)

}
