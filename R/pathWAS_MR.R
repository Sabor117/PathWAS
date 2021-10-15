#' Conduct Mendelian randomisation of clumped SNPs from pathway genes on end-point omics
#'
#' @description
#' Inputs aqnd combines aspects from your QTL SNPs, clumped SNPs, pathway genes, and omics summary
#' stats (plus flip status) and conducts MR between these, outputting the MR for model.
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
#' The mr_mvinput function creates an output specifically for use by the mr_mvlasso function which is then used
#' to run the MR of the QTL SNPs against the omics SNPs. Any of the exposures (I.e. genes) which then have an individual P-value of < 0.05 are considered to be significant
#' exposures within this model. The function will then return both the output from mr_mvlasso and the list of
#' significant exposures (this is used for downstream filtering).
#' Additionally, the option has been provided to save the output at multiple stages of the analysis so that they
#' can be examined individually. The MR lasso input and output can be saved, along with the genes which are
#' actually used to create the model. If you would like to save these you have to define a location where you
#' wish the output RDS file to be saved. You also have to input a name for the end-point and pathway (only used
#' for file-naming purposes).
#'
#' @param genelist list. List of genes extracted from database for your pathway which overlap with your QTL data.
#' @param clumped_snps data frame. The data frame of clumped SNPs output from qtl_clumpR.
#' @param qtl_sumstats character. File location for the QTLs used (either divided by gene or complete data).
#' @param geneCol character. Name or number of column containing the name of the gene for each SNP in the QTL data. Default is "gene". Ensure the format of the name is the same as your gene list.
#' @param omics_snps data frame. Summary stats of end-point omcis, now including FLIP column frm omics_MungeR
#' @param omics_SNPCol character. Name/number of column containing the SNP ID in the omics data frame.
#' @param omics_BetaCol character. Name/number of column containing the effect size in the omics data frame.
#' @param omics_SECol character. Name/number of column containing the standard error in the omics data frame.
#' @param save_MRInput logical. Should the function save the input for the MR? Default is FALSE
#' @param save_MRInLoc character. If save_MRInput == TRUE, then provide a directory to save the .rds file to.
#' @param save_MROutput logical. Should the function save the output from the MR (as well as output it from the function)?
#' @param save_MROutLoc logical. If save_MROutput == TRUE, then provide a directory to save the .rds file to.
#' @param save_MRExps logical. Should the function save the list of gene exposures from the MR (as well as output it from the function)?
#' @param save_MRExpsLoc logical. If save_MRExps == TRUE, then provide a directory to save the .rds file to.
#' @param end_point character. Name of end point protein (can be any string). Used for saving purposes if any save variable is set to TRUE.
#' @param path_select character. Name of pathway (can be any string). Used for saving purposes if any save variable is set to TRUE.
#'
#' @examples
#' ## Havimg created a list of genes, clumped SNPs and munged your Omics SNPs you can input them in the following way, along with adding the option to save the MR output:
#' pathWAS_MR(genelist, clumped_snps, qtl_sumstats = "/opt/localdir/gene_qtls/file1_sumstats_$$$.tsv", geneCol = "gene_ensembl", omics_snps, save_MROutput = TRUE, save_MROutLoc = "/opt/localdir/mr_outputs/", path_select = "nod_signalling", end_point = "IL18")
#'
#' @import data.table MendelianRandomization
#'
#' @export
pathWAS_MR = function(genelist,
                         clumped_snps,
                         qtl_sumstats,
                         geneCol = "gene",
                         omics_snps,
                         omics_SNPCol= "rsid", omics_BetaCol = "beta1", omics_SECol = "se",
                         save_MRInput = FALSE,
                         save_MRInLoc = NULL,
                         save_MROutput = FALSE,
                         save_MROutLoc = NULL,
                         save_MRExps = FALSE,
                         save_MRExpsLoc = NULL,
                         end_point = NULL, path_select = NULL,
                         verbose = TRUE
                        ) {

  if (verbose == TRUE){

    heading("Nobody expects the MR inquisition... Our chief weapon is surprise, surprise and difficulty explaining the concept.")

  }

  if (any(c(save_MRInput, save_MROutput, save_MRExps))){

    if (is.null(end_point) || is.null(path_select)){

      stop("Cannot save MR files without input for end-point and pathway.\n====\n")

    }

    selected_saves = which(c(save_MRInput, save_MROutput, save_MRExps))
    selected_saves_names = c("MR Input", "MR Output", "MR exposures")

    save_file_locs = c(is.null(save_MRInLoc), is.null(save_MROutLoc), is.null(save_MRExpsLoc))

    for (savecheck in selected_saves){

      if (save_file_locs[savecheck]){

        stop(paste0("Please enter a save file location for ",
                    selected_saves_names[savecheck],
                    ".\n====\n")
        )

      }
    }
  }

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

  if (verbose == TRUE){

    heading("MR input made. Running MR.")

  }

  if (save_MRInput == TRUE) {

    saveRDS(mr_input,
            file = paste0(save_MRInLoc, end_point, "_", path_select, "_MR_input.rds"))


  }

  mr_lasso_res = MendelianRandomization::mr_mvlasso(mr_input)

  mr_lasso_res@Exposure = colnames(snp_beta_matrix)

  if (save_MROutput == TRUE) {

    saveRDS(mr_lasso_res,
            file = paste0(save_MROutLoc, end_point, "_", path_select, "_MR_output.rds"))

  }

  if (verbose == TRUE){

    heading("MR complete.")

  }

  mr_lasso_names = list(genes = colnames(snp_beta_matrix), rsids = rownames(snp_beta_matrix))

  if (save_MRExps == TRUE) {

    saveRDS(mr_lasso_names,
            file = paste0(save_MRExpsLoc, end_point, "_", path_select, "_MR_exps.rds"))

  }

  mr_sig_index = mr_lasso_res$Pvalue < 0.05
  sig_mr_genelist = mr_lasso_res@Exposure[mr_sig_index]

  return(c(mr_lasso_res, sig_mr_genelist))

}
