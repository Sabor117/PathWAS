#' Conduct Elasticnet Mendelian randomisation of clumped SNPs from pathway genes on end-point omics
#'
#' @description
#' Variant of the pathWAS_MR function, performs elasticnet mendelian randomisation on input SNPs and PRS from PathWAS
#' pipeline. Outputs an elasticnet object which can be used for further downstream analysis.
#'
#' @details
#' This function requires the input of numerous pieces of data from throughout the PathWAS pipeline. The list
#' of genes (acquired from searching for the pathway genes) is needed to select specific QTL SNPs for conducting
#' the MR, this list must exclude the end-point gene. A data frame containing the SNPs produced by clumping (qtl_clumpR) is needed for pruning the QTLs by
#' clumped SNPs. A file location of the data frame of your QTL SNPs including a column with the named gene which must match the format
#' of the input gene list (I.e. HGNC format), this can either be the total available sumstats or can point to
#' multiple files divided up by gene. This should also include all SNPs for every gene and not just significant ones.
#' Lastly the data frame of your end-point omics. This must have first been through the omics_MungeR function
#' to add the "FLIP" column (you will now also be required to define columns for the effect size and standard error).
#' The QTL file location can either be a total file of all genes and SNPs in your data, or it can be individual
#' files divided up by gene (for easier/quicker reading). If it is the latter, then the file name must be in the
#' following format: "/opt/dir/eqtl_data_$$$.tsv" where "$$$" will be the name of one of the gene.
#' With all of this provided this function first creates two matrices of SNPs which overlap between the clumped
#' SNPs and the QTL SNPs. This will provide a matrix of SNPs against genes and provide the standard error and beta
#' of each. In the cases where there is no data for the clumped SNP for one of the genes, the standard error and
#' beta are both set to extremely low (1 and 0.0000001 respectively).
#' These matrices are used as the input for mr_mvinput from the MendelianRandomization R package (which is required
#' to run this function). This function also requires the betas and SEs from the omics SNPs, with the betas now
#' aligned to the QTLs using the newly created "FLIP" column.
#' You have the option of providing a previously created MR input object which allows most of this to be skipped
#' and the function will immediately create the elastic net output.
#'
#' @param genelist list. List of genes extracted from database for your pathway which overlap with your QTL data.
#' @param clumped_snps data frame. The data frame of clumped SNPs output from qtl_clumpR.
#' @param qtl_sumstats character. File location for the QTLs used (either divided by gene or complete data).
#' @param geneCol character. Name or number of column containing the name of the gene for each SNP in the QTL data. Default is "gene". Ensure the format of the name is the same as your gene list.
#' @param omics_snps data frame. Summary stats of end-point omcis, now including FLIP column frm omics_MungeR
#' @param omics_SNPCol character. Name/number of column containing the SNP ID in the omics data frame.
#' @param omics_BetaCol character. Name/number of column containing the effect size in the omics data frame.
#' @param omics_SECol character. Name/number of column containing the standard error in the omics data frame.
#' @param save_MRInput MR Input object. Contains a matrix of exposure SNPs to betas, exposure SNPs to SEs and then a list of outcome SNP betas and SEs. Default is NULL so function will create one for you.
#' @param end_point character. Name of end point protein (can be any string). Used for saving purposes if any save variable is set to TRUE.
#' @param path_select character. Name of pathway (can be any string). Used for saving purposes if any save variable is set to TRUE.
#' @param verbose logical. Enables chatty error and progress reporting.
#'
#' @import data.table MendelianRandomization
#'
#' @export
pathWAS_enet = function(genelist,
                      clumped_snps,
                      qtl_sumstats,
                      MRInput = NULL,
                      geneCol = "gene",
                      omics_snps,
                      omics_SNPCol= "rsid", omics_BetaCol = "beta1", omics_SECol = "se",
                      end_point = NULL, path_select = NULL,
                      verbose = TRUE
                      ) {

  if (verbose == TRUE){

    heading("African or European Mendelian Randomisation?")

  }

  if (is.null(MRInput)){

    clumped_snplist = clumped_snps$rsid
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

      matrix_col = colnames(snp_beta_matrix)

      snp_beta_matrix = as.matrix(data.frame(snp_beta_matrix[omics_snps[, omics_SNPCol],]))
      snp_se_matrix = as.matrix(data.frame(snp_se_matrix[omics_snps[, omics_SNPCol],]))

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

      snp_beta_matrix = snp_beta_matrix[omics_snps[, omics_SNPCol],]
      snp_se_matrix = snp_se_matrix[omics_snps[, omics_SNPCol],]

    }

    if (verbose == TRUE){

      heading("Matrices made. Creating elastic net MR input.")

    }

    omics_betas = omics_snps[, omics_BetaCol] * omics_snps$FLIP
    omics_se = omics_snps[, omics_SECol]

    mr_input = MendelianRandomization::mr_mvinput(bx = snp_beta_matrix,
                                                  bxse = snp_se_matrix,
                                                  by = omics_betas,
                                                  byse = omics_se)

    if (verbose == TRUE){

      heading("MR input made. Running MR.")

    }
  } else {

    mr_input = MRInput

    path_cohort_ovgenes = unique(genelist[!(genelist %in% end_point)])

  }

  mr_enet_res = glmnet_enet_mr(mr_input)

  if (verbose == TRUE){

    cat("\nElastic Net residuals\n___________\n\n")

    print(mr_enet_rest)

  }

  if (is.null(MRInput)){

    if (length(path_cohort_ovgenes) != length(mr_enet_res@Exposure)){

      warning("Input gene list has a different number of exposures than in MR input object.")

    }

    mr_enet_res@Exposure = colnames(snp_beta_matrix)

  } else {

    mr_enet_res@Exposure = path_cohort_ovgenes

  }

  if (verbose == TRUE){

    heading("Elastic net complete.")

  }

  return(mr_enet_res)

}
