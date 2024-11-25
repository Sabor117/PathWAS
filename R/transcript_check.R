#' Check whether a gene is transcribed via the pathway (I.e. via interaction with DNA)
#'
#' @description
#' transcript_check returns a single row of a data frame, containing the interactions leading to it in KEGG.
#'
#' @details
#' This function uses a defined KEGG KGML, pathway name, and Entrz gene. It searches the KGML for this gene
#' (only in nodes which the gene is at the end of an interaction) and then returns those rows to the function
#' user. If the "type" of interaction is GErel, this is defined by KEGG as "gene expression interaction,
#' indicating relation of transcription factor and target gene product" and therefore, this is a gene product
#' which is transcribed through the pathway.
#'
#' @param pathway character. KEGG pathway ID (e.g.: hsa05131).
#' @param gene_entrez character or numeric. Entrez gene ID.
#' @param kegg_file character. Full file location of a KEGG KGML file.
#'
#' @import KEGGgraph
#'
#' @export
transcript_check = function(pathway,
                             gene_entrez,
                             kegg_file = "./"){

  cat(paste0("Checking whether your end-point (hsa:", gene_entrez, ") is a gene expressed via pathway ", pathway, ".\n"))
  cat("Come check out the expression in the system! Help! Help! I'm being expressed!\n\n")

  pathway_frame = parseKGML2DataFrame(kegg_file, reactions = TRUE)

  gene_row = pathway_frame[pathway_frame$to == paste0("hsa:", gene_entrez),]

  print(gene_row)

  return(gene_row)
}
