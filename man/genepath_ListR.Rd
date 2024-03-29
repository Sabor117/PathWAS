% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genepath_ListR.R
\name{genepath_ListR}
\alias{genepath_ListR}
\title{Create pathway gene-lists for all GTEx tissues}
\usage{
genepath_ListR(
  gene,
  pathway,
  tissue = NULL,
  write = TRUE,
  genelistDir = getwd(),
  hsapien_mart,
  transcriptFile,
  kgmlDir = getwd(),
  delete_tmps = FALSE,
  cut_paths = -1,
  keep_routes = TRUE
)
}
\arguments{
\item{gene}{character or numeric (Entrez ID)}

\item{pathway}{character. KEGG pathway ID (path:hsa5200).}

\item{tissue}{character. Name of specific tissue for creating gene-list (must be as in GTEx TPMs file). Default is NULL which works for all tissues.}

\item{genelistDir}{directory. Directory to search for and save created data frames. Default is current working directory.}

\item{hsapien_mart}{data frame. Data frame of gene names including at least entrezgene_id, external_gene_name.}

\item{transcriptFile}{data frame. GTEx TPMs file read as data frame.}

\item{kgmlDir}{directory. Directory to search for and save KEGG KGML files. Default is current working directory.}

\item{delete_tmps}{logical. Optionally delete KGMLs after creating gene list.}

\item{cut_paths}{numeric. Cut-off for length of simple paths (default is -1 for all_smple_paths). Add threshold to reduce computing time.}
}
\description{
genepath_ListR returns a data frame of genes for a given pathway in each tissue available in GTEx TPMs.
}
\details{
This is a higher order function consisting of a run-through the creation of a gene-list for any given pathway
in KEGG. The function searches in the provided directory (genelistDir) for a file with the name in the format
<ENDPOINT>_<PATHWAY>_genelist.txt. If it finds one, it reads it in. If it does not find one it will create
one and then save it to the specified directory.

This defaults to create a gene-list for every tissue in GTEx (including a "Complete" gene-list which
incorporates every gene in the pathway with no filtering based on TPMs) for the chosen end-point + pathway
combination.

This function requires several inputs. It requires the GTEx TPM files (downloaded from https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz),
and it requires a BioMart file (from https://www.ensembl.org/biomart/martview/ or using package biomaRt)
which includes the following columns: entrezgene_id, external_gene_name.
}
