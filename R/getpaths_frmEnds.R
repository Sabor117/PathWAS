#' Get (KEGG) biological pathways from end point
#'
#' @description
#' getpaths_frmEnds returns a list of (KEGG) pathways from an input gene end-point.
#'
#' @details
#' This function requires the input of a single gene (in Entrez ID format) from this it will
#' search the KEGG database for any pathways which include the gene. This is done by converting the
#' Entrez ID into a KEGG human gene ID (E.g. 1234 becomes hsa:1234).
#' The pathways in which the gene is present are then filtered to only include those where the gene
#' is an "end-point" for the pathway: I.e. it has one or more inward-facing edges and no outward-facing
#' edges.
#' These pathways are those which are passed to return()
#' If no pathways are found or the gene is present in pathways but is not an end-point, the function
#' will inform the user and will return nothing.
#'
#' @param gene character or numeric (Entrez ID)
#'
#' @examples
#' ## Search for pathways in KEGG which have IL18 as an end-point.
#' getpaths_frmEnds(3606)
#'
#' @import KEGGREST KEGGgraph KEGGlincs splitstackshape
#'
#' @export
getpaths_frmEnds = function(gene_entrez){

  cat(paste0("Obtaining pathways for: ", gene_entrez, ".\n"))
  cat("It's only a wafer-thin mint, sir...\n\n")

  geneKEGG = paste0("hsa:", gene_entrez)

  if (length(gene_entrez) == 0) {

    heading(paste0("Gene ", gene_entrez, " not present in KEGG: ", geneKEGG))
    cat("\n\n")

    return()

  }

  pathfind = as.list(KEGGREST::keggLink("pathway", geneKEGG))

  print(head(pathfind))
  cat("\n\n")

  keeps = list()

  if (length(pathfind) > 0){

    for (endPI in 1:length(pathfind)){

      pathway_check = pathfind[[endPI]]

      tmp_fl = tempfile()

      cat(paste0("Lookup KEGG pathway for ", gene_entrez," + ", pathway_check, ".\n\n"))

      pathway_check = gsub("path:", "", pathway_check)

      kgml_file = KEGGlincs::get_KGML(pathway_check)

      kgml_mappings = KEGGlincs::expand_KEGG_mappings(kgml_file, FALSE)

      kgml_edges = KEGGlincs::expand_KEGG_edges(kgml_file, kgml_mappings)
      pathway_edges = KEGGlincs::edge_mapping_info(kgml_edges)

      cat("Lookup successful. Creating edges frame.\n\n")

      expanded_edges = data.frame(in_node = pathway_edges$entry1symbol,
                                  out_node = pathway_edges$entry2symbol)

      expanded_edges = cSplit(expanded_edges, "in_node", ",", "long")
      expanded_edges = cSplit(expanded_edges, "out_node", ",", "long")

      path_indegrees = nrow(expanded_edges[expanded_edges$out_node == gene_entrez,])
      path_outdegrees = nrow(expanded_edges[expanded_edges$in_node == gene_entrez,])

      if (path_outdegrees == 0){

        cat(paste0("Gene ", gene_entrez, " end-point for pathway: ", pathway_check, ".\n"))

        if (path_indegrees > 0){

          cat(paste0("With ", path_indegrees, " in-edges. Kept."))
          cat("\n\n===\n\n")

          keeps = c(keeps, pathway_check)

        } else {

          cat(paste0("With ", path_indegrees, " in-edges. Not kept."))
          cat("\n\n===\n\n")

        }

      } else {

        cat(paste0("Gene ", gene_entrez, " NOT end-point for pathway: ", pathway_check, "."))
        cat("\n\n===\n\n")
      }
    }
  }

  keeps = unlist(keeps)

  if (length(keeps) > 0){

    keeps = data.frame(pathway = keeps)
    row.names(keeps) = c()

    heading("getpaths_frmEnds complete. List of pathways created.")
    print(head(keeps))
    cat("\n\n\n=================\n\n")

    return(keeps)

  } else {

    heading("getpaths_frmEnds complete.")
    cat(paste0("Gene ", gene_entrez, " is not an end-point for any KEGG pathways."))
    cat("\n\n\n=================\n\n")

  }
}
