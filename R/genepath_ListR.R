#' Create pathway gene-lists for all GTEx tissues
#'
#' @description
#' genepath_ListR returns a data frame of genes for a given pathway in each tissue available in GTEx TPMs.
#'
#' @details
#' This is a higher order function consisting of a run-through the creation of a gene-list for any given pathway
#' in KEGG. The function searches in the provided directory (genelistDir) for a file with the name in the format
#' <ENDPOINT>_<PATHWAY>_genelist.txt. If it finds one, it reads it in. If it does not find one it will create
#' one and then save it to the specified directory.
#'
#' This defaults to create a gene-list for every tissue in GTEx (including a "Complete" gene-list which
#' incorporates every gene in the pathway with no filtering based on TPMs) for the chosen end-point + pathway
#' combination.
#'
#' This function requires several inputs. It requires the GTEx TPM files (downloaded from https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz),
#' and it requires a BioMart file (from https://www.ensembl.org/biomart/martview/ or using package biomaRt)
#' which includes the following columns: entrezgene_id, external_gene_name.
#'
#' @param gene character or numeric (Entrez ID)
#' @param pathway character. KEGG pathway ID (path:hsa5200).
#' @param tissue character. Name of specific tissue for creating gene-list (must be as in GTEx TPMs file). Default is NULL which works for all tissues.
#' @param genelistDir directory. Directory to search for and save created data frames. Default is current working directory.
#' @param hsapien_mart data frame. Data frame of gene names including at least entrezgene_id, external_gene_name.
#' @param transcriptFile data frame. GTEx TPMs file read as data frame.
#' @param kgmlDir directory. Directory to search for and save KEGG KGML files. Default is current working directory.
#' @param delete_tmps logical. Optionally delete KGMLs after creating gene list.
#' @param cut_paths numeric. Cut-off for length of simple paths (default is -1 for all_smple_paths). Add threshold to reduce computing time.
#'
#'
#' @import data.table
#'
#' @export
genepath_ListR = function(gene, pathway, tissue = NULL,
                          write = TRUE,
                          genelistDir = getwd(),
                          hsapien_mart,
                          transcriptFile,
                          kgmlDir = getwd(),
                          delete_tmps = FALSE,
                          cut_paths = -1,
                          keep_routes = TRUE
                          ) {

  readmart = fread(hsapien_mart,
                   data.table = FALSE)

  geneName = unique(readmart$external_gene_name[readmart$entrezgene_id %in% gene])
  pathName = gsub("path:", "", pathway)

  if (length(geneName) > 1){

    warning(paste0("genepath_ListR WARN1: Multiple instances of gene name (", geneName, ") found in BiomaRt. Attempting to prune by HGNC symbol."))

    cat(paste0("\n\nNames for gene found:\n\n"))
    print(geneName)

    refinedMart = readmart[readmart$entrezgene_id %in% gene,]
    refinedMart = refinedMart[!(refinedMart$hgnc_symbol == ""),]

    geneName = unique(refinedMart$external_gene_name[refinedMart$entrezgene_id %in% gene])

    if(length(geneName) == 1){

      cat(paste0("\n\nGene name selected: ", geneName,"\n==========\n\n"))

    } else{

      stop(paste0("genepath_ListR Error1: Multiple instances of gene name found in BiomaRt. Unable to prune by HGNC. Please check: ", geneName, "\n-----\n\n"))

    }
  }

  if (file.exists(paste0(genelistDir, geneName, "_", pathName, "_genelist.txt")) == FALSE){

    ### If no existing combination file is found, then it will create one:

    heading(paste0("No existing gene list/s for pathway + gene combo (", pathName, " + ", geneName, "). Making you one now."))

    ### Obtain list of genes from KEGG for pathway

    pathway_genes = kegg_geneglob(pathway, hsapien_mart)

    ### Obtain list of all simple pathways from KEGG pathway

    all_simple_paths = smple_paths(pathway, gene,
                                   hsapien_mart = readmart,
                                   saveDir = kgmlDir,
                                   delete_tmp = delete_tmps,
                                   cut_paths = cut_paths,
                                   keep_routes = keep_routes)

    if (length(all_simple_paths) == 0){

      heading("No simple paths to end-gene.")

      return()
    }

    all_smple_paths_genes = unique(unlist(all_simple_paths))

    ### Obtain genes specific to tissues within pathway
    ### Requires expression data (TPMs) from GTEx for each tissue
    ### tissue variable must be in GTEx format
    ### If variable is default null or "all" then run on all tissues

    gtex_tpms = data.table::fread(transcriptFile,
                        data.table = FALSE)

    ### If a tissue was specified then only create data frame for that tissue
    ### However, do not save the file (only save overall combinations in a file)

    if (is.null(tissue) || toupper(tissue) == "ALL"){

      heading("All tissues specified, creating a list of all tissues + gene-lists.")

      tissue_genes_list = get_tissuegenes(gtex_tpms, pathway_genes, all_simple_paths) ### tissue_genes_list outputs as a list of lists

      if (length(tissue_genes_list) == 0){

        heading("No genes-paths in any tissue. Only using Complete.")

        tissue_genes_frame = data.frame(tissue = "Complete", gene = all_smple_paths_genes)

        if (write == TRUE){

          heading("Write called as true, saving as data frame.")

          data.table::fwrite(tissue_genes_frame,
                             paste0(genelistDir, geneName, "_", pathName, "_genelist.txt"),
                             sep = "\t", quote = FALSE, row.names = FALSE)

        }

      } else {

        tissue_genes_frame = c()

        ### Loop to convert list into a data frame of two columns

        for (ncomb in names(tissue_genes_list)){

          tissue_genes_frame = rbind(tissue_genes_frame, cbind(ncomb, tissue_genes_list[[ncomb]]))

        }

        ### Changes column names to tissue + gene (I.e. every gene + tissue combination)

        colnames(tissue_genes_frame) = c("tissue", "gene")

        ### Add a set of genes for the complete pathway

        comp_path_frame = data.frame(tissue = "Complete", gene = all_smple_paths_genes)

        tissue_genes_frame = rbind(tissue_genes_frame, comp_path_frame)

        ### As this loop involves ALL combinations, it will save the data frame to a file
        ### EDIT: specifies directory

        if (write == TRUE){

          heading("Write called as true, saving as data frame.")

          data.table::fwrite(tissue_genes_frame,
                             paste0(genelistDir, geneName, "_", pathName, "_genelist.txt"),
                             sep = "\t", quote = FALSE, row.names = FALSE)

        }

      }
    } else { ### If no file exists but a tissue is specified:

      heading(paste0(tissue, " tissue specified, creating a gene list for it."))

      ### Creates tissue_gene list but only for specified tissue

      tissue_genes_list = get_tissuegenes(gtex_tpms, pathway_genes, all_simple_paths, tissue)

      if (length(tissue_genes_list) == 0){

        heading("No genes in selected tissue.")

      } else {

        ### Converts list into data frame

        for (ncomb in names(tissue_genes_list)){

          tissue_genes_frame = rbind(tissue_genes_frame, cbind(ncomb, tissue_genes_list[[ncomb]]))
        }

        ### Changes column names to tissue + gene (I.e. every gene + tissue combination)
        ### Does NOT save

        colnames(tissue_genes_frame) = c("tissue", "gene")

        heading("Results not saved to disk.")

      }
    }

    ### If a file already exists for this particular pathway

  } else if (file.exists(paste0(genelistDir, geneName, "_", pathName, "_genelist.txt")) == TRUE){

    heading("Pathway + end-point gene list exists. Reading it now.")

    ### Instead of creating list - simply reads in data frame of combinations

    all_tissue_genes = data.table::fread(paste0(genelistDir, geneName, "_", pathName, "_genelist.txt"),
                             data.table = FALSE)

    ### selects specific parts of the data frame based on whether a single tissue was specified

    if (is.null(tissue) || toupper(tissue) == "ALL"){

      heading("Creating gene lists for all tissues.")

      tissue_genes_frame = all_tissue_genes

    } else {

      heading("Creating gene lists for selected tissue.")

      tissue_genes_frame = all_tissue_genes[all_tissue_genes$tissue == tissue,]

    }
  }

  ### Making sure tissue_genes is in data.frame

  if (!(exists("tissue_genes_frame"))){

    heading("No possible gene lists for pathway + end-point combo.")
    return(NULL)

  } else {

    tissue_genes_frame = data.frame(tissue_genes_frame)

  }

  return(tissue_genes_frame)

}
