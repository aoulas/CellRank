# CellRank
## Update
Jun 27, 2023 (Version 1.0.0)

CellRank is now available via devtools installation. 

<div style='text-align: justify;'>
  Some very long text here.
  The text must span multiple lines in order to demonstrate that this works.mation to guide the choice of cell types from a scRNA-seq analysis that yield the most biologically meaningful results. Prior knowledge is incorporated in a standardized, structured manner, whereby a checklist is attained by querying MalaCards human disease database with a disease of interest. The checklist is comprised of pathways and drugs and optionally drug mode of actions (MOAs), associated with the disease. The user is prompted to “edit” this checklist by removing or adding terms (in the form of keywords) from the list of predefined terms. The user may also define de-novo, a set of keywords that best suit a hypothesis the user is interested in investigating (hypothesis-driven approach). Once the checklist is finalized a “mapping” step is performed. This is done initially against pathway enrichment results attained from analysing the scRNA-seq data. In addition, the user-
  bobbyhadz.com
</div>

## Capabilities
CellRank utilizes prior knowledge in combination with expert-user information to guide the choice of cell types from a scRNA-seq analysis that yield the most biologically meaningful results. Prior knowledge is incorporated in a standardized, structured manner, whereby a checklist is attained by querying MalaCards human disease database with a disease of interest. The checklist is comprised of pathways and drugs and optionally drug mode of actions (MOAs), associated with the disease. The user is prompted to “edit” this checklist by removing or adding terms (in the form of keywords) from the list of predefined terms. The user may also define de-novo, a set of keywords that best suit a hypothesis the user is interested in investigating (hypothesis-driven approach). Once the checklist is finalized a “mapping” step is performed. This is done initially against pathway enrichment results attained from analysing the scRNA-seq data. In addition, the user-selected checklist of drug names is mapped against DR result derived from the data. The analysis of scRNA-seq data has the capability to obtain pathways and repurposed drugs by comparing disease-control conditions for every cell type in the analysis. CellRank then uses the user defined information to highlight the cell types that generate the results that best “map” to the predefined prior knowledge provided by the expert user. The methodology is fully automated and a ranking is generated for all cell types in the analysis allowing the user to pinpoint the specific cells that are most prominently affected by the disease under study in accordance to the provided prior knowledge.  The output of our methodology provides an automated validation of the result and further makes it easier for researchers to interpret their findings. In addition, the results provide greater credence to de novo information as it also backed-up by prior knowledge for the disease under study. Novel information is obtained in the form of predicted pathways that are enriched for each cell type, as well as previously unreported repurposed drugs that are obtained from the differentially expressed genes (DEGs) between disease-control conditions.   

## Installation
CellRank R package can be easily installed from Github using devtools:

devtools::install_github("aoulas/CellRank")

Please make sure you have installed all the dependencies. See instruction below.

## Installation of  dependencies
### CRAN packages
Seurat, dplyr, patchwork, multtest, metap, ggplot2, cowplot, enrichR, gridExtra, ggpubr, RColorBrewer, crank, riverplot, rvest, stringr.
### Bioconductor packages
KEGGREST, GO.db, rWikiPathways, ReactomeContentService4R, msigdb.
### GitHub packages
SeuratDisk, CellChat.

Some users might have issues when installing CellChat pacakge due to different operating systems and new R version. Please check the following solutions:

Installation on [Windows](https://github.com/sqjin/CellChat/issues/5)
