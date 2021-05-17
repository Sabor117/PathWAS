#' Obtain all genes from a KEGG pathway
#'
#' @description
#' Returns all genes associated with input pathway (excluding directionality)
#'
#' @details
#' This function takes a pathway (possibly defined by the getpaths_frmEnds function) from KEGG and then finds
#' all genes associated with this pathway and outputs that.
#' Depending on whether or not you have an available BioMart reference table, this can be either as a list of Entrez IDs or as a data frame
#' of KEGG IDs (which is simply "hsa:ENTREZ-ID") and gene names.
#' A BioMart table can be obtained either by downloading a full dataset from https://www.ensembl.org/biomart/martview/ or through usage of
#' the biomaRt package. This table MUST have the entrezgene_id and external_gene_name columns (default column names in BioMart for Entrez
#' and gene names).
#'
#' @param pathway character. KEGG pathway ID.
#' @param in_mart data frame. A data frame of multiple columns which, at a minimum includes entrezgene_id and external_gene_name.
kegg_geneglob = function(pathway,
                         in_mart = NULL){

  cat(paste0("Obtaining KEGG gene list for: ", pathway, ".\n"))
  cat("Nobody expects the KEGG inquisition.\n\n")

  genelist = KEGGREST::keggLink(pathway)[,2] ### specifically selects second column (genes)
  genelist = genelist[grep("hsa", genelist)] ### specifies human genes from second column
  genelist = gsub("hsa:", "", genelist) ### convert KEGG IDs to Entrez

  if (!(is.null(in_mart))){

    ### Convert Entrez IDs in genelist into gene names
    ### WARNING: currently unable to deal with multiple Entrez IDs for different genes

    cat("Lookup Biomart.\n\n")

    id_map = data.table::fread(in_mart,
                   data.table = FALSE)

    ### EDIT: Column number depends on id_map used

    genelist_frame = data.frame(unique(id_map[id_map$entrezgene_id %in% genelist, c("entrezgene_id", "external_gene_name")]))

    ### Keep both KEGG ID and HGNC ID

    genelist_frame$entrezgene_id = paste0("hsa:", genelist_frame$entrezgene_id)

    colnames(genelist_frame) = c("kegg_id", "external_gene_name")

    genelist = genelist_frame

    heading("kegg_geneglob complete. Genelist created.")
    print(head(genelist))
    cat("\n\n")

    return(genelist)

  } else if (is.null(in_mart)){

    heading("kegg_geneglob complete. Genelist created.")
    print(head(genelist))
    cat("\n\n")

    return(genelist)

  }
}
