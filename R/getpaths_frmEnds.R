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
#' @import KEGGREST KEGGgraph KEGGlincs splitstackshape igraph
#'
#' @export
getpaths_frmEnds = function(gene_entrez,
                            hsapien_mart,
                            kgmlDir = getwd(),
                            saveKGML = TRUE){

  cat(paste0("Obtaining pathways for: ", gene_entrez, ".\n"))
  cat("It's only a wafer-thin mint, sir...\n\n")

  geneKEGG = paste0("hsa:", gene_entrez)

  if (length(gene_entrez) == 0) {

    heading(paste0("Gene ", gene_entrez, " not present in KEGG: ", geneKEGG))
    cat("\n\n")

    return()

  }

  gene_name = unique(hsapien_mart$external_gene_name[hsapien_mart$entrezgene_id %in% gene_entrez])

  if (length(gene_name) > 1){

    warning("genepath_ListR WARN1: Multiple instances of gene name found in BiomaRt.")

  }

  pathfind = as.list(KEGGREST::keggLink("pathway", geneKEGG))

  print(head(pathfind))
  cat("\n\n")

  keeps = list()

  if (length(pathfind) > 0){

    for (endPI in 1:length(pathfind)){

      pathway_check = pathfind[[endPI]]

      cat(paste0("Lookup KEGG pathway for ", gene_entrez," + ", pathway_check, ".\n\n"))

      pathway_check = gsub("path:", "", pathway_check)

      path_save_file = paste0(kgmlDir, pathway_check, "_kegg_file.kgml")

      tmp_fl = tempfile() ### Necessary for retrieveKGML

      pathway_kgml = try(KEGGgraph::retrieveKGML(pathway_check, ### Search for given pathway
                                                 organism = "hsa", ### Organism
                                                 destfile = tmp_fl, ### This is necessary for some reason
                                                 method = "wget", ### Utilises wget method
                                                 quiet = TRUE))

      if (file.exists(paste0(path_save_file))){

        cat(paste0("KGML file for pathway exists here: ", path_save_file, "\n\n"))

      } else {

        cat(paste0("Downloading KGML file for pathway here: ", path_save_file, "\n\n"))

        system(paste0("wget ", pathway_kgml, " -O ", path_save_file))

      }

      pathway_info = KEGGgraph::parseKGML2Graph(path_save_file, ### pathway kgml file
                                                expandGenes = TRUE, ### expand paralogue nodes
                                                genesOnly = FALSE) ### include connections to things which aren't genes

      pathway_table = igraph::as_long_data_frame(igraph::igraph.from.graphNEL(pathway_info))

      path_indegrees = nrow(pathway_table[pathway_table$to_name == geneKEGG,])
      path_outdegrees = nrow(pathway_table[pathway_table$from_name == geneKEGG,])

      if (!(path_outdegrees == 0)){

        cat(paste0("Gene ", gene_entrez, " NOT end-point for pathway: ", pathway_check, "."))
        cat("\n\n============\n\n")

      } else {

        cat(paste0("Gene ", gene_entrez, " end-point for pathway: ", pathway_check, ".\n"))

        if (path_indegrees == 0){

          cat(paste0("With ", path_indegrees, " in-edges. Not kept."))
          cat("\n\n============\n\n")

        } else if (path_indegrees > 0){

          cat(paste0("With ", path_indegrees, " in-edges."))
          cat("\n\n")
          cat("Looking up KEGGlincs.\n\n")

          kgml_file = try(KEGGlincs::get_KGML(pathway_check))

          if (!(class(kgml_file) == "KEGGPathway")){

            warning(paste0("getpaths_frmEnds WARN2: No KEGGlincs download error for: ", pathway_check, " with end-point: ",gene_entrez))

            next

          }

          if (all(is.na(kgml_file@nodes))){

            warning(paste0("getpaths_frmEnds WARN1: No KEGGlincs KGML nodes or edges for pathway: ", pathway_check))

            next

          }

          kgml_mappings = KEGGlincs::expand_KEGG_mappings(kgml_file, FALSE)

          complex_frame = kgml_mappings[grepl(gene_name, kgml_mappings$LABEL),]
          complex_status = any(grepl("Complex", complex_frame$LABEL))

          if (complex_status == FALSE){

            cat(paste0("Not a complex. VALID END POINT. Kept.\n\n============\n\n"))

            keeps = c(keeps, pathway_check)

          } else {

            complex_genes = complex_frame$LABEL[grepl("Complex", complex_frame$LABEL)][[1]]

            cat(paste0("End-point found as ", complex_genes, ". Confirming validity.\n\n"))

            complex_genes = unlist(strsplit(gsub(" Complex", "", complex_genes), ":"))

            complex_entrez_ids = unique(hsapien_mart$entrezgene_id[hsapien_mart$external_gene_name %in% complex_genes])

            complex_check = c()

            for (num_complex in 1:length(complex_entrez_ids)){

              currKEGG = paste0("hsa:", complex_entrez_ids[num_complex])

              path_outdegrees = nrow(pathway_table[pathway_table$from_name == currKEGG,])

              if (path_outdegrees > 0){

                complex_check = c(complex_check, TRUE)

              } else {

                complex_check = c(complex_check, FALSE)

              }
            }

            if (any(complex_check)){

              cat(paste0("Complex has external edges. NOT VALID END POINT.\n\n============\n\n"))

            } else {

              cat(paste0("Complex genes have no external edges. VALID END POINT. Kept\n\n============\n\n"))

              keeps = c(keeps, pathway_check)

            }
          }
        }
      }

      if (saveKGML == FALSE){

        system(paste0("rm ", path_save_file))

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
