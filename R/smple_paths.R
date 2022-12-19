#' Get all simple paths from KEGG pathway
#'
#' @description
#' smple_paths returns a list of all simple (one-route) pathways through a larger/more complex bioogical pathway to a chosen end point.
#'
#' @details
#' This function takes a pathway (possibly defined by the getpaths_frmEnds function) from KEGG and then finds
#' all simple "routes" that lead towards a chosen end-point. Therefore it searches for every route through the pathway
#' and searches for every gene included in those routes that lead towards the end-point.
#' Therefore this outputs only the genes from the pathway which are associated with the end-point, and excludes genes which
#' theoretically have no influence over the end-point.
#'
#' @param pathway character. KEGG pathway ID (e.g.: path:hsa05131).
#' @param gene_entrez character or numeric. Entrez gene ID.
#' @param keep_routes logical. Should the function output all simple pathways, or output a list of genes.
#' @param saveDir character. Directory for saving the KGML file downloaded from KEGG
#' @param delete_tmps logical. After running, delete the KGML downloaded or not.
#' @param hsapien_mart data frame. Data frame of gene names including at least entrezgene_id, external_gene_name.
#'
#' @examples
#' ## Find all simple routes in pathway hsa05131 (Shigellosis) leading to end point 3606 (IL18).
#' ## Output only the genes involved in those routes (I.e. do not keep the routes themselves).
#'
#' smple_paths(pathway = path:hsa05131, gene_entrez = 3606, keep_routes = FALSE)
#'
#' @import KEGGgraph igraph
#'
#' @export
smple_paths = function(pathway,
                       gene_entrez,
                       keep_routes = TRUE,
                       hsapien_mart,
                       saveDir = "./",
                       delete_tmps = FALSE){

  cat(paste0("Obtaining all simple paths list for: ", pathway, ".\n"))
  cat("Come Sergei!\n\n")

  path_check = pathway
  geneKEGG = paste0("hsa:", gene_entrez)
  geneNAME = unique(hsapien_mart$external_gene_name[hsapien_mart$entrezgene_id %in% gene_entrez])

  path_save_name = gsub("path:", "", path_check)
  path_save_file = paste0(saveDir, path_save_name, "_kegg_file.kgml")

  tmp_fl = tempfile() ### Necessary for retrieveKGML

  cat("Lookup KEGG.\n\n")

  pathway_kgml = try(KEGGgraph::retrieveKGML(path_check, ### Search for given pathway
                                  organism = "hsa", ### Organism
                                  destfile = tmp_fl, ### This is necessary for some reason
                                  method = "wget", ### Utilises wget method
                                  quiet = TRUE))


  if (file.exists(paste0(path_save_file))){

    cat(paste0("KGML file for pathway exists here: ", path_save_file, "\n\n"))

  } else {

    cat(paste0("Downloading KGML file for pathway here: ", path_save_file, "\n\n"))

    system(paste0("wget ", pathway_kgml, " -O ", saveDir, path_save_name, "_kegg_file.kgml"))

  }

  pathway_info = KEGGgraph::parseKGML2Graph(path_save_file, ### pathway kgml file
                                            expandGenes = TRUE, ### expand paralogue nodes
                                            genesOnly = FALSE) ### include connections to things which aren't genes

  pathway_table = igraph::as_long_data_frame(igraph::igraph.from.graphNEL(pathway_info))

  cat("Download and read successful.\n\n")

  last_layer = pathway_table[pathway_table$to_name == geneKEGG,]

  simplified_pathway_table = last_layer

  connected_nodes = unique(last_layer$from_name)
  gene_check = connected_nodes

  repetitions = 0

  repeat {

    repetitions = repetitions + 1

    workingNode = connected_nodes[repetitions]

    new_layer = pathway_table[pathway_table$to_name == workingNode,]
    simplified_pathway_table = unique(rbind(simplified_pathway_table, new_layer))

    new_nodes = unique(new_layer$from_name)
    connected_nodes = unique(c(connected_nodes, new_nodes))

    if (length(connected_nodes) <= repetitions){

      break

    }
  }

  ### Make simplified igraph object from simple table
  ### Convert names back to KEGG names

  connected_nodes = c(connected_nodes, geneKEGG)

  simplified_igraph = igraph::graph_from_data_frame(simplified_pathway_table,
                                                    directed = TRUE,
                                                    vertices = NULL)

  vertice_match = unique(data.frame(node = c(simplified_pathway_table$from, simplified_pathway_table$to),
                                    name = c(simplified_pathway_table$from_name, simplified_pathway_table$to_name)))

  igraph::vertex_attr(simplified_igraph)$name = vertice_match$name[match(igraph::vertex_attr(simplified_igraph)$name, as.character(vertice_match$node))]

  ### Define "starting" genes by those which are never in the "to" column

  start_genes = which(!(simplified_pathway_table$from %in% simplified_pathway_table$to))

  if (length(start_genes) == 0){

    warning("Error smple_paths1: Pathway contains no 'starting' genes. Simplified pathway may be a loop.")

    return(connected_nodes)

  }

  dir_paths = list() ### Empty list

  ### For every starter gene in pathway:

  for (nstart in 1:length(start_genes)){

    ### Select start gene

    nstart_gene = simplified_pathway_table$from_name[start_genes[nstart]]

    ### all_simple_paths function from igraph
    ### Finds/outputs all straight line connections from selected start gene
    ### To defined end point

    smple_path_n = igraph::all_simple_paths(simplified_igraph,
                                            nstart_gene,
                                            to = which(igraph::vertex_attr(simplified_igraph)$name == geneKEGG)) ### Select vertice of end point

    ### Creating list of all simple paths

    if (length(smple_path_n) != 0){ ### If simple path does not exist from start point, skip back to start

      dir_paths_part = list() ### Initialise part of list

      for (nsmples in 1:length(smple_path_n)){

        ### For every item in the "part" list, split into individual IDs and make into a normal list
        ### Normal "part" list then added to overall list

        curr_smple = igraph::as_ids(smple_path_n[[nsmples]])

        dir_paths_part[[nsmples]] = curr_smple

      }

      ### Combine into one list

      dir_paths = unique(c(dir_paths, dir_paths_part))

    }
  }

  dir_paths = unique(dir_paths)

  heading("smple_paths complete. All simple paths created. Simples.")
  print(head(dir_paths))
  cat("\n\n")

  if (delete_tmps == TRUE){

    cat("Deleting downloaded KGML file.\n\n")

    system(paste0("rm ", path_save_file))

  }

  if (keep_routes == TRUE){

    return(dir_paths)

  } else {

    return(unique(unlist(dir_paths)))

  }
}
