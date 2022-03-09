# PathWAS

:arrow\_right: ASHG 2021 poster is available
[here](doc/).

## Contents


## Contributors

The majority of the code for PathWAS was written by Sebastian May-Wilson while working with Dr. Nicola Pirastu.

The package relies on numerous other packages including: MendelianRandomization, KEGGgraph, KEGGrest and ieugwasr. It also includes code written by Dr. Verena Zuber.


## Overview

*PathWAS* is an R package designed to attempt to predict pathway functionality from multi-omics data. It attempts to create a polygenic risk score (PRS) for a given pathway (extracted from online databases) using a measured -omics (usually proteomics) end-point as a proxy for total functionality. This model can then be applied to PheWAS methods to try and find complex traits associated with your pathway/s of interest.

This is a multi-step process which ideally can be run from top to bottom with the correct data in the right format.

This README will describe the rationale of the package, then list the main and important functions of PathWAS, followed finally by a demonstration of how to run the package from top to bottom.


## Rationale

The concept of *PathWAS* is that if we assume a biological pathway is a series of chemical reactions relying on the expression of different genes/proteins, then increased expression of these genes will lead to increased (or decreased) activity at that step of the pathway. Therefore in it's simplest terms the combined activity of the pathway can be estimated by the cumulative expression of each gene in the pathway.

Genes relating to a pathway can be either extracted from an online database (e.g. KEGG, Reactome, etc) or personally curated.

We obtain the expression levels of the pathway genes from QTLs. Again, these can be either in-house QTLs or from numerous online sources. The PathWAS package was created and tested using both eQTLs from GTEx v.8 and eQTLgen. In theory, however, it may be even more useful to use pQTLs if a dataset is available with enough gene coverage (an important aspect of PathWAS is including as many genes from a pathway as possible and at present most pQTL datasets are limited to smaller numbers of proteins than genes in eQTL sets).

These QTLs are then used in a multi-variable MR against your selected end-point -omics. PathWAS was created and tested using in-house proteomics from the SCALLOP consortium and from the ORCADES cohort, but theoretically if you have the GWAS and measurements of another form of omics (e.g. metabolomics or even transcriptomics) you can use these. In this step the SNPs from your QTLs are clumped and then used to create an MR input object. These SNPs are used as instrumental variables for the MR (with the genes as exposures) against the SNPs from your selected omics.


## Pre-requisites of PathWAS

This section contains a list of files and inputs (and the format of the various files, where necessary) for the complete running of PathWAS. Some files are highlighted as being optional.

- **-Omics end-point GWAS** - A GWAS of one or more genes selected as a potential end-point for PathWAS. The omics in question can be metabolomics/transcriptomics but both have drawbacks over usage of proteomics. The GWAS must contain the following information: rsID, beta/effect, standard error of effect, effect allele (A1), other allele (A0). This data must be in data frame format, column names are unimportant. 


## Primary functions and requirements of PathWAS

- **`getpaths_frmEnds()`** - KEGG-SPECIFIC FUNCTION. Takes an input gene or protein (in Entrez gene ID format) and searches the KEGG database for pathways containing the gene. Then outputs all of those pathways for which the protein is an end-point.

- **`smple_paths()`** - KEGG-SPECIFIC FUNCTION. Takes an input pathway and an Entrez end-point gene and refines the selected pathway to a list of simple pathways connected to the end-point. Collectively these simple pathways represent the part of the pathway that is relevant to the chosen end-point, and a list of genes can be extracted from it.

- **`get_tissuegenes()`** - Using the list of simple paths created by `smple_paths()` this function utilizes GTEx TPM expression data to refine the selected genes and simple pathways further to create tissue-specific gene lists. The requires a local downloaded copy of the GTEx TPMs file (found here: https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz)

- **`genepath_ListR()`** - This is a wrapper for the `smple_paths()` and `get_tissuegenes()` functions. Inputting an end-point gene (Entrez ID), a pathway name (KEGG pathway ID) and then the GTEx TPMs file + a BioMart file containing gene translations



