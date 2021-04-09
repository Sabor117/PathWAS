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
#' @param gene character (Entrez ID)
#'
#' @example
#' ## Search for pathways in KEGG which have IL18 as an end-point.
#' getpaths_frmEnds(3606)
getpaths_frmEnds = function(gene_entrez){

  cat(paste0("Obtaining pathways for: ", gene_entrez, ".\n"))
  cat("It's only a wafer-thin mint, sir...\n\n")

  require(KEGGREST)
  require(KEGGgraph)

  cat("Lookup successful.\n\n")

  geneKEGG = paste0("hsa:", gene_entrez)

  if (length(gene_entrez) == 0) {

    heading(paste0("Gene ", gene_entrez, " not present in KEGG: ", geneKEGG))
    cat("\n\n")

    return()

  }

  pathfind = as.list(keggLink("pathway", geneKEGG))

  print(head(pathfind))
  cat("\n\n")

  keeps = list()

  if (length(pathfind) > 0){

    for (endPI in 1:length(pathfind)){

      pathway_check = pathfind[[endPI]]

      tmp_fl = tempfile()

      cat(paste0("Lookup KEGG pathway for ", gene_entrez," + ", pathway_check, ".\n\n"))

      pathway_kgml =	try(retrieveKGML(pathway_check,
                                      organism = "hsa",
                                      destfile = tmp_fl,
                                      method = "wget",
                                      quiet = TRUE))

      pathway_info = parseKGML2Graph(pathway_kgml,
                                     expandGenes = TRUE,
                                     genesOnly = FALSE)

      cat("Lookup successful.\n\n")

      path_edges = KEGGgraph::edges(pathway_info)
      path_outdegrees = sapply(KEGGgraph::edges(pathway_info), length)

      path_indegrees = sapply(KEGGgraph::inEdges(pathway_info), length)

      if (path_outdegrees[[geneKEGG]] == 0){

        cat(paste0("Gene ", gene_entrez, " end-point for pathway: ", pathway_check, ".\n"))

        if (path_indegrees[[geneKEGG]] > 0){

          cat(paste0("With ", path_indegrees[[geneKEGG]], " in-edges. Kept."))
          cat("\n\n")

          keeps = c(keeps, pathway_check)

        } else {

          cat(paste0("With ", path_indegrees[[geneKEGG]], " in-edges. Not kept."))
          cat("\n\n")

        }

      } else {

        cat(paste0("Gene ", gene_entrez, " NOT end-point for pathway: ", pathway_check, "."))
        cat("\n\n")
      }
    }
  }

  if (length(keeps) > 0){

    keeps = t(as.data.frame(keeps))
    rownames(keeps) = c()
    colnames(keeps) = c("pathway")

    heading("getpaths_frmEnds complete. List of pathways created.")
    print(head(keeps))
    cat("\n\n")

    return(keeps)

  } else {

    heading("getpaths_frmEnds complete.")
    cat(paste0("Gene ", gene_entrez, " is not an end-point for any KEGG pathways."))
    cat("\n\n\n")

  }
}
