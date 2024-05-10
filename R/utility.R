

#' Euclidean distances.
#'
#' Compute Euclidean distances between two points. X and Y coordinates
#' must be given for each point. 
#'
#' @param coord1 A 2 element numeric vector with coordinates for point 1. 
#' @param coord2 A 2 element numeric vector with coordinates for point 2. 
#'
#' @return The euclidean distance between the two given coordinates.
#' 
#' @export 
#'
eucDist = function(coord1, coord2) {
  # Setup.
  res = 0
  for (i in seq_len(length(coord1))) {
    res = res + ((coord1[i] - coord2[i]) ^ 2)
    eucDist = sqrt(res)
  }
  return(eucDist)
}


#' Find cluster centroids.
#'
#' Compute the cluster centroids for a given set of cell type labels. They could
#' be cell type annotation or any cluster identification number. 
#'
#' @param seuratObj A Seurat object with dim reduced embeddings.
#' @param embeddingKey A character indicating the name of the embedding (e.g., pca).
#' @param groupLabels A character indicating the label by which to group the
#'  cells for computing the centroids. 
#' @param dims A numeric vector with the dimensions of the embedding across which we
#'  wish to compute the cluster centroids. Default is to use only the the first 2 dims. 
#'  
#' @return A dataframe with centroids computed across each PC or embedding 
#'  column vector. 
#'  
#' @export
#'
findClusterCentroids = function(
    seuratObj, 
    embeddingKey, 
    groupLabels,
    dims = c(1, 2)
    ) {
  # Check provided embedding is present.
  if (!embeddingKey %in% names(seuratObj@reductions)) {
    msg = paste0(embeddingKey, " not present in Seurat object embeddings!")
    stop(msg)
  }
  embeddings = Embeddings(seuratObj, reduction = embeddingKey)
  centroidData = aggregate(
    embeddings, 
    by = list(seuratObj[[]][, groupLabels]),
    FUN = mean
  )
  centroidData = tibble::column_to_rownames(
    centroidData, 
    var = "Group.1"
    )
  # Subset to desired dimensions. 
  centroidData = centroidData[, dims]
  return(centroidData)
}


# Distance between centroids. 
# calcCentroidDist = function(centroids, label) {
#   targetClusterCoords = centroids[label, ]
#   distCentroids = apply(
#     centroids,
#     MARGIN = 1,
#     FUN = function(x) eucDist(targetClusterCoords, x) 
#     )
#   return(distCentroids)
# }


#' Compute centroid distances. 
#'
#' Find centroids and compute distance of a target annotation centroid 
#' to centroids of other clusters.
#' 
#' @param seuratobj A Seurat object.
#' @param embeddingKey A character indicating the name of the embedding used
#'  for comuting centroids. 
#' @param groupLabels A character indicating the name of the cell type annotation
#'  across which to compute centroids. 
#' @param targetLabel the target annotation that we want to comapre other
#'  annotations with. Euclidean distance of target annotation to itself should
#'  be zero of course. 
#' 
#' @return A dataframe with euclidean distance of a target cell annotation centroid
#'  to the centroids of other clusters or annotations. 
#' @export
#'
calcCentroidDist = function(
    seuratObj, 
    embeddingKey, 
    groupLabels, 
    targetLabel
    ) {
  # Find cluster centroids.
  centroids = findClusterCentroids(
    seuratObj = seuratObj,
    embeddingKey = embeddingKey,
    groupLabels = groupLabels
    )
  # Compute of target annotation to other cluster centroids. 
  targetClusterCoords = centroids[targetLabel, ]
  centroidDists = apply(
    centroids,
    MARGIN = 1,
    FUN = function(x) eucDist(targetClusterCoords, x) 
  )
  # Prep df for plotting. 
  centroidDistsDf = data.table::rbindlist(centroidDists)
  colnames(centroidDistsDf)[1] = "euclidean_distance"
  centroidDistsDf$targetAnno = targetLabel
  centroidDistsDf$cluster = names(centroidDists)
    
  return(centroidDistsDf)
}












