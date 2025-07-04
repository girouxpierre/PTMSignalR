# PTMSignalR (BulkSignalR extension for PTM integration)

## Overview
This extension of BulkSignalR ([here](https://github.com/jcolinge/BulkSignalR)) is used to integrate post-translational modifications (PTMs) in ligand-receptor (L-R) interactions inference. It can work with bulk or single-cell expression (transcriptomics/proteomics) data and can integrate phosphorylation, ubiquitination and glycosylation data. Other modifications databases can be added by the user. Potential L-R interactions are taken from the LRdb database, which is included in our other package SingleCellSignalR, available from Bioconductor [here](https://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html).
Inferences rely on a statistical model linking potential L-R interactions with biological pathways from Reactome or biological processes from GO.
A number of visualization and data summary functions are proposed to help navigating the predicted interactions.

## Installation
``
# BulkSignalR and this extension are not included in BioConductor yet.
# Installation goes via GitHub:
# install.packages("devtools")
devtools::install_github("girouxpierre/PTMSignalR",build_vignettes = TRUE)

# To read the vignette
# browseVignettes("PTMSignalR")
``
