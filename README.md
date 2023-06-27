# CellRank
Update
JUN 27, 2023 (Version 1.0.0)

CellRank is now available vai devtools installation. 

Capabilities
CellRank utilize prior knowledge in combination with expert-user information to guide the choice of cell types from a scRNA-seq analysis that yield the most biologically meaningful results. Prior knowledge is incorporated in a standardized, structured manner, whereby a checklist is attained by querying MalaCards human disease database (6) with a disease of interest (e.g., Covid-19). The checklist is comprised of pathways and drugs and optionally drug mode of actions (MOAs), associated with the disease. The user is prompted to “edit” this checklist by removing or adding terms (in   the form of keywords) from the list of predefined terms. The user may also define de-novo, a set of keywords that best suit a hypothesis the user is interested in investigating (hypothesis-driven approach). Once the checklist is finalized a “mapping” step is performed. This is done initially against pathway enrichment results attained from analysing the scRNA-seq data. In addition, the user-selected checklist of drug names is mapped against DR result derived from the data. The analysis of scRNA-seq data has the capability to obtain pathways and repurposed drugs by comparing disease-control conditions for every cell type in the analysis. Our methodology then uses the user defined information to highlight the cell types that generate the results that best “map” to the predefined prior knowledge provided by the expert user. The methodology is fully automated and a ranking is generated for all cell types in the analysis allowing the user to pinpoint the specific cells that are most prominently affected by the disease under study in accordance to the provided prior knowledge (see Figure 1).  The output of our methodology provides an automated validation of the result and further makes it easier for researchers to interpret their findings. In addition, the results provide greater credence to de novo information as it also backed-up by prior knowledge for the disease under study. Novel information is obtained in the form of predicted pathways that are enriched for each cell type, as well as previously unreported repurposed drugs that are obtained from the differentially expressed genes (DEGs) between disease-control conditions.   

Installation
CellRank R package can be easily installed from Github using devtools:

devtools::install_github("aoulas/CellRank")
Please make sure you have installed all the dependencies. See instruction below.

Installation of  dependencies
dplyr,Seurat,SeuratDisk,patchwork,multtest,metap,ggplot2,cowplot,enrichR,gridExtra,ggpubr,RColorBrewer,crank,KEGGREST,GO.db,rWikiPathways,ReactomeContentService4R,msigdb,riverplot,CellChat,dplyr,rvest,stringr

Some users might have issues when installing CellChat pacakge due to different operating systems and new R version. Please check the following solutions:

Installation on Windows
