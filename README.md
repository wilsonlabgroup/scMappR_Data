
# single-cell mapper (scMappR) data repository

### Dustin Sokolowski: dustin-dot-sokolowski-at-sickkids-dot-ca

### Date: 09/28/2020


## Description


### Primary utility of the scMappR_Data R package.
This repository complements the primary scMappR package by providing all of the data files required to run the different analyses in scMappR. Every function within the main package has a `rda_path` parameter. This parameter allows users to set where the scMappR_Data package is downloaded in their local machine. If the `scMappR_Data` pacakge is not downloaded and a file is required, then scMappR will temporarily download and load that file into R using the `downloader` R package.

### Primary function of scMappR

*Description of scMappR R package adapted from pre-print*. The single cell mapper (scMappR) R package contains a suite of bioinformatic tools that provide experimentally relevant cell-type specific information to a list of differentially expressed genes (DEG). 
The function `scMappR_and_pathway_analysis` reranks DEGs to generate cell-type specificity scores called cell-weighted fold-changes (cwFold-change). cwFold-changes are the original fold-changes identified with bulk differential analysis that are then scaled by the cell-type specificity of the DEG and the cell-type proportion of the inputted samples.

#### "PBMC_scMappR_and_pathway_analysis_example.rda"

The "PBMC_scMappR_and_pathway_analysis_example.rda" file loads the PBMC_example object of class "list". This file contains the sex differences, signature matrix, and normalized bulk RNA-seq matrix from the Monaco et al., 2019 study (Figure 3 of the main figure).
* bulk_DE_cors: Gene, adjusted P-value, and Log2FC
* bulk_normalized: CPM normalized counts.
* odds_ratio_in: signature matrix of different immune cell-types.

####  "Preoptic_region_example.rda"
The "Preoptic_region_example.rda" file loads the POA_example of class "list". This file contains the cell-type markers, p-value signature matrix, and odds-ratio signature matrix of a Pre-optic area scRNA-seq dataset (Moffitt et al., 2018).
* POA_generes: A list of cell-type markers as an output of the `FindMarkers` function in `Seurat`.
* POA_OR_signature: The signature matrix of the fold-changes generated from the `2^(avg_logFC)*sign(ave_logFC)` of cell-type markers.
* POA_OR_signature: The signature matrix of the ranks generated from the -1*log10(p_val_adj)*sign(p_val_adj) of cell-type markers.
#### "Signature_matrices_OR.rda"
The "Signature_matrices_OR.rda" file loads the OddsRatioSignature of class "list". This file contains the fold-change signature matrices of all pre-processed scRNA-seq datasets.
#### "Signature_matrices_Pval.rda"
The "Signature_matrices_Pval.rda" file loads the RankValueSignature of class "list". This file contains the fold-change signature matrices of all pre-processed scRNA-seq datasets.
#### "Signature_matrices_cellLabel.rda"
The "Signature_matrices_cellLabel.rda" file loads the celltypeLabels of class "list". This file contains the top (up-to) 30 cell-type markers and estimated cell-type labels from the `CellMarker` and `Panglao` datasets using the Fisher's exact test and GSVA methods.
#### "TopGenes.rda"
The "TopGenes.rda" file loads the TopGenesList of class "list". This file contains the cell-type markers from each study, tissue, and cell-type.
#### "bioMart_ortholog_human_mouse.rda"
The "bioMart_ortholog_human_mouse.rda" file loads the bioMart_orthologs of class "data.frame". This file contains the 1-1 orthologs between mouse and human gene symbols.
#### "cell_preferences_categorized.rda
The "cell_preferences_categorized.rda" file loads the cell_preference_final of class "list". This file contains 14 tissue categories. If a researcher chooses one of the available cell-type labelling methods (e.g. brain). Then scMappR will pick labels from this category.
#### "cell_process_example.rda"
The "cell_process_example.rda" file loads the sm of class "dgCMatrix". This file contains a matrix of 22k genes and 500 cells from scRNA-seq retinal data. This is to work with the `process_dgTMatrix_lists`.
#### "do_valle_duraes_2020_rnaseq.rda"
The "do_valle_duraes_2020_rnaseq.rda" file loads the do_valle_duraes_2020_rnaseq of class "list". 
* counts: raw counts of bulk kidney RNA-seq from do Valle Duraes 2020.
* phenotype: data frame of the sample information from the bulk kidney RNA-seq samples.
* Signature: a gene by cell-type matrix of the fold-change that a gene is in each cell-type. The signature matrix is from the Tabula Muris 2018 paper.
#### "human_cell_markers"
The "human_cell_markers" file loads the gmt_list object that will label cell-types (human gene symbols) from different sources.
* gmt_both: cell-type marker database from the `CellMarker` and `PanglaoDB databases`.
* gmt_cellmarker: cell-type marker database from the `CellMarker` database.
* gmt_panglao: cell-type marker database from the `Panglao` database.
* gmt_gobp: list of human gene ontologies.
* gmt_subtype: cell-type marker database of neurons specifically.
#### "mouse_cell_markers"
The "mouse_cell_markers" file loads the gmt_list object that will label cell-types (human gene symbols) from different sources.
* gmt_both: cell-type marker database from the `CellMarker` and `PanglaoDB databases`.
* gmt_cellmarker: cell-type marker database from the `CellMarker` database.
* gmt_panglao: cell-type marker database from the `Panglao` database.
* gmt_gobp: list of human gene ontologies.
* gmt_subtype: cell-type marker database of neurons specifically.
#### "human_tissue_celltype_scMappR.rda"
The "human_tissue_celltype_scMappR.rda": cell-type marker database of human gene symbols listing human cell-types.
#### "mouse_tissue_celltype_scMappR.rda"
The "mouse_tissue_celltype_scMappR.rda": cell-type marker database of mouse gene symbols listing mouse cell-types.




