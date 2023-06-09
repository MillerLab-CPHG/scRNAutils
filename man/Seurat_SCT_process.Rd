% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_analysis.R
\name{Seurat_SCT_process}
\alias{Seurat_SCT_process}
\title{Downstream processing of counts matrix}
\usage{
Seurat_SCT_process(
  seurat_obj,
  seurat_filter = FALSE,
  min_UMI = 200,
  max_UMI = 20000,
  min_features = 200,
  max_features = 4000,
  max_mt_percent = 10,
  max_hb_percent = 5,
  library_id,
  study_id,
  arterial_origin = NULL,
  disease_status = NULL,
  sex = NULL
)
}
\arguments{
\item{seurat_obj}{An input seurat object.}

\item{seurat_filter}{A boolean indicating whether do we wish to filter. 
Set to FALSE in pre-processing round for doublet removal.}

\item{min_UMI}{A numeric indicating min number of UMIs to keep.}

\item{max_UMI}{A numeric indicating max number of UMIs to keep.}

\item{min_features}{A numeric indicating the minimum number of genes cells should express to be included.}

\item{max_features}{A numeric indicating the maximum number of genes cells should express to be included.}

\item{max_mt_percent}{A numeric indicating the max percentage of reads mapped to the mito genome.}

\item{max_hb_percent}{A numeric indicating the max percentage of reads mapped to hemoglobin genes.}

\item{library_id}{A character indicating name of the sample.}

\item{arterial_origin}{A character indicating arterial bed of the library.}

\item{disease_status}{A character indicating the disease status of the library (e.g., non-lesion, lesion).}

\item{sex}{A character string indicating sex of the subject.}

\item{study_name}{A character indicating name of the study.}
}
\value{
A seurat object that has been SCT normalized, nearest neighbors graph and UMAP embeddings.
}
\description{
Internal helper function to standardize filtering, SCT normalization
and dimensionality reduction. This function will also add
some important metadata variables.
}
