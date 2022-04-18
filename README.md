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

This section contains a list of files and inputs (and the format of the various files, where necessary) for the complete running of PathWAS. Some files/steps are highlighted as being optional.

- **-Omics end-point GWAS** - A GWAS of one or more genes selected as a potential end-point for PathWAS. The omics in question can be metabolomics/transcriptomics but both have drawbacks over usage of proteomics. The GWAS must contain the following information: rsID, beta/effect, standard error of effect, effect allele (A1), other allele (A0). This data must be in data frame format, column names are unimportant. 

- **cis-QTLs for genes** - Summary stats for cis-QTLs which will be used for the genes of a pathway. These can be either pQTLs or eQTLs however we recommend using data which covers as many genes as possible so that the pathway scores created will incorporate as many of the genes involved as possible. Large possible data sets include GTEx and eQTLgen.

- **Local plink installation** - This is required for clumping the QTL SNPs.

- **Local BioMart data frame of genes and name translations** - Important for numerous steps within PathWAS and given that many different sources of data will use different nomenclature. It may not be required, but often is.

- **List of pathway genes - OPTIONAL** - If you are examining a pre-determined or self-curated pathway you may skip the initial steps of PathWAS (which are designed to search for pathways from a given end-point). If you already have a list of genes making up your pathway leading in to your selected end-point -omics then you can skip these functions.

- **GTEx tissue TPMs - OPTIONAL** - Used by a specific function which will refine your gene list by expression within tissues based on the GTEx TPMs file. This is an optional step.

- **Second, independent measurement of the same -omics** - Used specifically in the final prediction stage of PathWAS. This is not a technical requirement for running PathWAS, but is one we recommend if testing multiple pathways.



## Primary functions and requirements of PathWAS

This does list not include functions used within other PathWAS functions.

- **`getpaths_frmEnds()`** - KEGG-SPECIFIC FUNCTION. Takes an input gene or protein (in Entrez gene ID format) and searches the KEGG database for pathways containing the gene. Then outputs all of those pathways for which the protein is an end-point. This is the first step required for running PathWAS.

- **`smple_paths()`** - KEGG-SPECIFIC FUNCTION. Takes an input pathway and an Entrez end-point gene and refines the selected pathway to a list of simple pathways connected to the end-point. Collectively these simple pathways represent the part of the pathway that is relevant to the chosen end-point, and a list of genes can be extracted from it.

- **`get_tissuegenes()`** - Using the list of simple paths created by `smple_paths()` this function utilizes GTEx TPM expression data to refine the selected genes and simple pathways further to create tissue-specific gene lists. The requires a local downloaded copy of the GTEx TPMs file (v8 found here: https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz).

- **`genepath_ListR()`** - This is a wrapper for the `smple_paths()` and `get_tissuegenes()` functions. Inputting an end-point gene (Entrez ID), a pathway name (KEGG pathway ID) and then the GTEx TPMs file + a BioMart file containing gene translations. This provides the `get_tissuegenes()` output (a data frame of tissue and each gene in the pathway within that tissue, plus a "complete" pathway list). This is the second step required for running PathWAS.

- **`qtl_clumpR()`** - A clumping function. This requires the input of your genes (this can be of the whole pathway output of `genepath_ListR()` or it can be a tissue-specific list) and your QTL SNPs from summary stats (this can be either all SNPs or only significant SNPs, however we recommend using significant associations). This function is a wrapper for the ieugwasr clumping function `ld_clump()` and as such requires a local installation of plink which is used for clumping the SNPs based on LD. This will output a data frame (in the same format as your input summary stats) of the clumped SNPs for your pathway genes. This is the third step required for running PathWAS.

- **`omics_MungeR()`** - Requires the data frame of clumped SNPs from `qtl_clumpR()` and the summary stats for the SNPs of the GWAS from your end-point -omics. This step will examine both data frames and tests whether the betas of the QTL SNPs need to be flipped, based on the alleles. It will output the omics SNPs data frame (refined to only the clumped SNPs) with an additional "flip" column. This is the fourth step required for running PathWAS.

- **`pathWAS_MR()`** - This step requires most pieces of data created thus far in PathWAS. It require the lit of gene names, the clumped SNPs, the omics SNPs (with the added "flip" column) and the end-point protein name. This step will then run the MR of the clumped SNPs against the omics SNPs to provide the gene effect on the pathway. This step also allows the option of saving the various MR inputs and outputs (for repeated use of the function). This is the fifth (and technically final) step of PathWAS. From this you will have a list of genes and their effects on the pathway based on an end-point. These MR effects can then be used to create PRS for each gene to create the overall pathway score.

- **`pathWAS_enet()`** - This step is functionally very similar to `pathWAS_MR()` but also applies an elastic net penalised regression (function written by Dr. Verena Zuber <v.zuber(at)imperial.ac.uk>). This is an additional and optional version of the PathWAS MR which did not yield results as robust as the normal MR.

- **`pathWAS_predictR()`** - This is a technically optional step of PathWAS, but one which we recommend carrying out and which was used in the PathWAS paper methodology. This requires an additional, prediction set of -omics data for the same end-point used to create the MR scores. You will need to use your initial set of -omics to create PRS for the end-point and then also have the measurements of the same -omic in the second data set. Inputting both of these, along with your MR effect sizes and the function will test how well your gene scores created from your first set of -omics predict the values of the second set. This step is not required for PathWAS (as the intention is not to have pathway PRS which can predict the end-point protein) but is a reasonable method of testing which pathways from your data are the best.








