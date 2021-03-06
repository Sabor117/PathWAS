#' Obtain a tissue-specific list of genes for a pathway
#'
#' @description
#' Returns a series of gene lists for either all tissues in the GTEx transcriptomics dataset or for a specific tissue.
#'
#' @details
#' This is a higher level function for creating one or multiple gene-lists for tissue-specific expression of genes. This is based on
#' the GTEx Transcripts Per Million data (which can be found here:
#' https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz).
#' From this data we define a threshold (default = 1) constituting "above basal expression" of the gene. Any gene above this threshold
#' is considered to be "expressed" in the defined tissue. The function then uses the linear paths provided from the smple_paths function and
#' excludes any simple path which has at one or more genes with an expression level below the defined threshold (under the assumption
#' that with one of the links in the chain missing, this is no longer a viable route for pathway functionality). Thus the output from this
#' is not just only the genes which are expressed in a specific tissue, but also theoretically will exclude those genes which can only
#' influence the desired end-point through genes which are defined as not being expressed within that tissue.
#'
#' This function requires several inputs. It requires the GTEx TPM files, the input gene-list must come in the format of a column with
#' kegg_id (this can be updated to be ENTREZ IDs) and external_gene_name (from BioMart). It also requires the output from smple_paths()
#' under the assumption of having kept the connections between each gene and not just the names of the genes.
#'
#' There default running option for this function provides the gene list for all tissues in the GTEx TPM files, however individual tissues
#' can be selected. There is currently no option for selecting multiple or a list of tissues (this may be updated down the line).
#'
#' @param exprn_data data frame. The GTEx TPM file having been read in as a data frame.
#' @param in_genelist data frame. A table of all genes in the pathway including a column with kegg_id (hsa:<ENTREZ>) and external_gene_name (gene name) as defined by BioMart
#' @param all_simple_paths list. The list of all "simple" linear paths through a pathway to a chosen end-point.
#' @param in_tissue character. Either NULL (and output a gene list for all tissues in GTEx) or a single specific tissue (the name must match the file names for the tissue in GTEx).
#' @param tpm_threshold numeric. Threshold defining whether a gene is expressed or not. Default = 1
#'
#' @import data.table
#'
#' @export
get_tissuegenes = function(exprn_data,
                           in_genelist,
                           all_simple_paths,
                           in_tissue = NULL,
                           tpm_threshold = 1){

  cat(paste0("Obtaining tissue-specific gene lists.\n"))
  cat("The Black Knight ALWAYS triumphs.\n\n")

  ### Subset all GTEx TPMs by pathway genelist

  pathway_genes_tissues = exprn_data[exprn_data$Description %in% in_genelist$external_gene_name,]

  ### In order to change row names to gene names - necessary to remove duplicate gene names

  duplicate_rmv = data.table::data.table(pathway_genes_tissues) ### convert to data table

  duplicate_rmv$sums = duplicate_rmv[, apply(.SD, 1, sum), .SDcols = sapply(duplicate_rmv, is.numeric)] ### creates a sum column of all tissue transcript values
  duplicate_rmv = duplicate_rmv[order(abs(sums), decreasing = T)][!duplicated(Description)] ### Entire table ordered by sums and then duplicates are removed if they are lower on the list
  duplicate_rmv$sums = NULL ### sums column removed

  pathway_genes_tissues = as.data.frame(duplicate_rmv)

  rownames(pathway_genes_tissues) = pathway_genes_tissues$Description ### Change rownames to gene names
  pathway_genes_tissues[,1:2] = NULL ### Remove gene name + ID columns

  ### Convert data.frame to matrix

  tissue_present = as.matrix(pathway_genes_tissues)
  tissue_present = matrix(as.numeric(tissue_present >= tpm_threshold), ### converts all values less than 1 into 0
                          nrow = nrow(pathway_genes_tissues), ### based on number of rows of original dataframe
                          byrow = F) ### dunno why

  ### Changes row and column names

  colnames(tissue_present) = colnames(pathway_genes_tissues)
  rownames(tissue_present) = rownames(pathway_genes_tissues)

  ### Create an empty matrix of size y by x
  ### Where y = number of genes in pathway
  ### And x = number of simple pathways
  ### Rownames changed to KEGG IDs

  simple_paths_matrix = matrix(nrow = nrow(in_genelist),
                               ncol = length(all_simple_paths))
  rownames(simple_paths_matrix) = in_genelist$kegg_id

  ### For each simple pathway

  cat("Creating simple paths matrix\n")

  for (nsimples in 1:length(all_simple_paths)){

    ### Single simple path selected
    ### Every gene present in pathway is ticked as "1" in previously empty matrix

    tester_path = all_simple_paths[[nsimples]]
    tester_path_col = as.numeric(rownames(simple_paths_matrix) %in% tester_path)
    simple_paths_matrix[,nsimples] = tester_path_col

  }

  cat("Matrix built.\n")

  ### Column names set to null (I.e. 1, 2, 3... n) - for simple paths

  colnames(simple_paths_matrix) = NULL

  ### Matrix is transposed for matrix multipliction

  simple_paths_matrix = t(as.matrix(simple_paths_matrix))

  ### KEGG IDs changed to HGNC symbols to match other matrix
  ### Based on genelist

  collist = in_genelist$external_gene_name[which(colnames(simple_paths_matrix) == in_genelist$kegg_id)]
  colnames(simple_paths_matrix) = collist

  ### Column order (I.e. genes) changed to match order of other matrix

  colorder = rownames(tissue_present)
  simple_paths_matrix = simple_paths_matrix[,colorder]

  ### Matrix multiplication

  tissue_path_matrix = simple_paths_matrix %*% tissue_present

  ### Number of genes in each simple path

  count_smpaths = lengths(all_simple_paths)

  ### Copy of tissue * path matrix with gene count in a column

  tissue_path_matrix_c = data.frame(tissue_path_matrix, count_smpaths)

  for (ncount in 1:nrow(tissue_path_matrix_c)){

    ### For every value in dataframe (matrix) compare with corresponding value in count column
    ### If value matches (I.e. tissue has every gene from simple path expressed)
    ### Change value in matrix to TRUE (1)

    tissue_path_matrix_c[ncount, -ncol(tissue_path_matrix_c)] = ifelse(tissue_path_matrix_c[ncount, 1:ncol(tissue_path_matrix_c)-1] == tissue_path_matrix_c[ncount, ncol(tissue_path_matrix_c)],
                                                                       TRUE,
                                                                       FALSE)
  }

  ### Convert dataframe back into matrix
  ### Removes gene count column

  tissue_path_matrix2 = as.matrix(tissue_path_matrix_c[,-ncol(tissue_path_matrix_c)])

  ### Initialise output list

  alltissue_genes = list()

  for (ntissue in 1:ncol(tissue_path_matrix2)){

    ### For every tissue (I.e. column)
    ### Select column specifically
    ### Find tissue name

    tiss_col = tissue_path_matrix2[,ntissue]
    tissue_head = colnames(tissue_path_matrix2)[ntissue]

    ### Initalise minilist

    tissue_genelist = c()

    for (numpath in 1:length(tiss_col)){

      if (tiss_col[numpath] == 1){

        ### For every simple path, compare with tissue * path matrix
        ### If value of corresponding cell = 1 (I.e. simple path has all genes expressed)
        ### Add genes from simple path to genelist

        tissue_genelist = c(tissue_genelist, all_simple_paths[[numpath]])

      }
    }

    ### Remove duplicate genes from genelist

    tissue_genelist = unique(tissue_genelist)

    ### If genelist contains no genes (I.e. non expressed in tissue)\
    ### Leave as 0
    ### Otherwise create list of tissues + all genes

    if (is.null(tissue_genelist) == FALSE){

      alltissue_genes[[ntissue]] = tissue_genelist

    } else {

      alltissue_genes[[ntissue]] = 0

    }

    ### Change names of list to match tissue names

    names(alltissue_genes)[ntissue] = tissue_head

  }
  ### Remove tissues for which there is no genes

  alltissue_genes = alltissue_genes[lengths(alltissue_genes) > 1]

  if (is.null(in_tissue)){

    ### Output all tissue + gene combinations

    heading("Tissue genelists created.")
    print(head(alltissue_genes))
    cat("\n\n")

    return(alltissue_genes)

  } else {

    ### Output specific tissue + gene combinations

    heading(paste0("Tissue genelist created for ", in_tissue, "."))
    print(head(alltissue_genes[in_tissue]))
    cat("\n\n")

    return(alltissue_genes[in_tissue])
  }
}
