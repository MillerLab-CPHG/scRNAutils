# scRNAutils

`scRNAutils` is a package that provides a comprehensive pipeline for processing single-cell RNA-seq data. Many of the functions are designed to complement gold-standard normalization and dimensionality reduction functionality within Seurat. Some features from the package include:

-   Robust QC including removal of doublets and ambient RNA with `scDblFinder` and `decontX`.

-   Adaptive thresholding QC metrics along with Median Absolute Deviation (MAD) outlier detection.

-   Automated selection of PCs explaining a given amount of variance (Default is 90%).

-   Several QC plots and summary statistics, including the N most highly expressed genes in the data and library complexity plots.

-   Optimization of clustering resolution using silhouette coefficients.

-   Other utility functions such as computing the correlation of a target gene across the entire transcriptome or distance of cell type cluster centroids.

### Installation

The package can be installed from Github with the following:

```         
if (!requireNamespace("remotes")) install.packages("remotes")
library(remotes)
remotes::install_github(MillerLab-CPHG/scRNAutils)
```

For more detailed documentation in how to run the main analysis pipeline as well as the expected inputs/outputs, please see the tutorial inside the `vignettes` folder.

**Cite scRNAutils:**

[Mosquera JV, Auguste G, Wong D, Turner AW, Hodonsky CJ, Alvarez-Yela AC, Song Y, Cheng Q, Lino Cardenas CL, Theofilatos K, Bos M, Kavousi M, Peyser PA, Mayr M, Kovacic JC, Björkegren JLM, Malhotra R, Stukenberg PT, Finn AV, van der Laan SW, Zang C, Sheffield NC, Miller CL. Integrative single-cell meta-analysis reveals disease-relevant vascular cell states and markers in human atherosclerosis. Cell Rep. 2023 Nov 28;42(11):113380. doi: 10.1016/j.celrep.2023.113380. Epub 2023 Nov 10. PMID: 37950869.](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)01392-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS221112472301392X%3Fshowall%3Dtrue)
