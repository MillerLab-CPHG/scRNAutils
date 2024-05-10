# scRNAutils
`scRNAutils` is a package that provides comprehensive pipeline for processing single-cell RNA-seq data. Many of the functions are designed to complement gold-standard normalization and dimensionality reduction functionality within Seurat. Besides adaptive and robust QC, the package optimizes clustering resolution parameters as well as selection of the N PCs explaining a given amount of variance in the data in the main analysis pipeline. We also offer some handy utility functions such as computing distances of a target cell annotation centroid to other annotation centroids. 

### Installation

The package can be installed from Github with:
`devtools::install_github(MillerLab-CPHG/scRNAutils)`

For more detailed documentation, please see the tutorial inside the `vignettes` folder.

**Cite scRNAutils:**
Mosquera JV, Auguste G, Wong D, Turner AW, Hodonsky CJ, Alvarez-Yela AC, Song Y, Cheng Q, Lino Cardenas CL, Theofilatos K, Bos M, Kavousi M, Peyser PA, Mayr M, Kovacic JC, Björkegren JLM, Malhotra R, Stukenberg PT, Finn AV, van der Laan SW, Zang C, Sheffield NC, Miller CL. Integrative single-cell meta-analysis reveals disease-relevant vascular cell states and markers in human atherosclerosis. Cell Rep. 2023 Nov 28;42(11):113380. doi: 10.1016/j.celrep.2023.113380. Epub 2023 Nov 10. PMID: 37950869.
