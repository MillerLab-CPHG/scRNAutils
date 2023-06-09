% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_expression_on_pseudotime}
\alias{plot_expression_on_pseudotime}
\title{Model gene expression across pseudotime.}
\usage{
plot_expression_on_pseudotime(
  cds_ordered,
  cell_annotations,
  genes,
  facet_wrap_plot = FALSE,
  pseudotime_boundary = 30
)
}
\arguments{
\item{genes}{A Vector with the genes to be plotted across pseudotime.}

\item{facet_wrap_plot}{A boolean indicating whether all genes shoud be shown in one or multiple plots.}

\item{pseudotime_boundary}{A numeric indicating the upper boundary of pseudotime to be plotted.}

\item{cds_ordered:}{A monocle3 cds object with cells ordered across pseudotime.}

\item{cell_annotations:}{A vector with cell types of interest for plotting gene expression.}
}
\value{
A ggplot object with one or multiple panels showing gene expression across pseudotime values.
}
\description{
This function will plot normalized expression values across pseudotime
}
