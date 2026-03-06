# PTMSignalR (BulkSignalR extension for PTM integration)

## Overview
This extension of BulkSignalR ([here](https://github.com/jcolinge/BulkSignalR)) is used to integrate post-translational modifications (PTMs) in ligand-receptor (L-R) interactions inference. It can work with bulk or single-cell expression (transcriptomics/proteomics) data and can integrate phosphorylation, ubiquitination and glycosylation data. Other modifications databases can be added by the user. Potential L-R interactions are taken from the LRdb database, which is included in our other package SingleCellSignalR, available from Bioconductor [here](https://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html).
Inferences rely on a statistical model linking potential L-R interactions with biological pathways from Reactome or biological processes from GO.
A number of visualization and data summary functions are proposed to help navigating the predicted interactions.

Please cite: Post-Transcriptional Modification Integration for Ligand–Receptor Cellular Network Inference. Giroux, et al. Mol Cell Proteomics, 2026.

## Installation
```
# BulkSignalR and this extension are not included in BioConductor yet.
# Installation goes via GitHub:
# install.packages("devtools")
devtools::install_github("girouxpierre/PTMSignalR",build_vignettes = TRUE)

# To read the vignette
# browseVignettes("PTMSignalR")
```

## Edit March 2026
Two different correction methods were implement to calculate final P-values in a more rigorous way. The original score obtained by the three P-values multiplication (L.R x R.T x R.PTM) can be replaced by corrected P-value obtained using Fisher’s method for P-value aggregation or using the method from Breitwieser et al. calculating the null-distribution of the product of independent biological sample P-values (Breitwieser et al., J Proteome Res, 2011).
The user can set the parameter "correction" when using the inference function:
```
# By default no correction is applied ("none") but the user can select one of the two correction methods.
initialInference(bsrdm.comp, "comp.name", correction = c("fisher", "breitwieser", "none"))

```
