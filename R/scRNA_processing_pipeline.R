

#' Process scRNA libraries 
#' 
#' A function that will handle QC (doublet and ambient RNA removal) and processing for individual scRNA libraries
#' 
#' @param count_matrix sparse count matrix required as initial input. As of now, make it accept 
#' a filtered_feature_bc_matrix with a Seurat::Read10X or a txt file. Maybe this is better
#' handled out of the function since we'd also be dealing with h5 files. Maybe we can do this later
#' but not priority right now
#' @param libraryID A character indicating the id for the library.
#' @param studyID A character string indicating the study that generated the data.
#' @param minRes A numeric indicating lower boundary for sil analysis
#' @param maxRes A numeric indicating upper boundary for sil analysis
#' @param dblFindIter A numeric indicating how many iterations we need to run for getting consensus doublets
#' @param queryGenes A character vector with genes of interest for checking expression before and after QC.
#' @param minUMI A numeric indicating min number of UMIs to keep.
#' @param maxUMI A numeric indicating max number of UMIs to keep.
#' @param minFeatures A numeric indicating the minimum number of genes cells should express to be included.
#' @param maxFeatures A numeric indicating the maximum number of genes cells should express to be included.
#' @param maxMtPercent A numeric indicating the max percentage of reads mapped to the mito genome. 
#' @param maxHbPercent A numeric indicating the max percentage of reads mapped to hemoglobin genes. 
#' @param arterialOrigin A character string indicating vascualr bed of the data.
#' @param diseaseStatus A character string indicating disease status of the library (e.g., lesion vs non-lesion)
#' @param sex A character string indicating sex of the subject. 
#' @return A list where the first object is the processed Seurat obj and remaining objects are UMAP plots before and after QC. 
#' @export
#' @examples 
#' seurat_outs = doitall(counts_matrix, "rpe004", "pan_et_al")

doItAll = function(
  countsMatrix, 
  libraryID,
  studyID,
  minRes = 0.3, 
  maxRes = 1.9, 
  dblFindIter = 3,
  setAutoThreshold = FALSE,
  regressMitoRibo = FALSE,
  rmMitoVarGenes = FALSE,
  rmRiboVarGenes = FALSE,
  clusteringAlg = "louvain",
  queryGenes = NULL,
  minUMI = 500,
  maxUMI = 20000,
  minFeatures = 200,
  maxFeatures = 4000,
  maxMtPercent = 10,
  maxHbPercent = 5,
  arterialOrigin = NULL,
  diseaseStatus = NULL,
  sex = NULL
  ){
  
  # Validate counts matrix input
  if (!methods::is(countsMatrix, "dgCMatrix")) {
    if(methods::is(countsMatrix, "Matrix")) {
      countsMatrix = as.sparse(countsMatrix)
    } else {
      msg = paste("Your counts matrix should be provided in matrix or 
                  sparse matrix format (dgCMatrix)!. Got ", 
                  class(countsMatrix))
      stop(msg)
    }
  }
  
  # How many rows and columns do this matrix have?
  message(paste("Your prefiltered count matrix has", 
                dim(countsMatrix)[1], "rows and",
                dim(countsMatrix)[2], "columns"))
  
  # QC (removal of doublets)
  # Prelim round of clustering
  seuratObj = CreateSeuratObject(counts = countsMatrix,
                                 project = libraryID,
                                 min.cells = 10,
                                 min.features = 200)
  
  # Do a first round of clustering for doublet removal
  # Set the seurat_filter arg to FALSE. We dont need to remove low quality cells
  # for this first round of clustering 
  seuratObj = seuratSCTprocess(seuratObj = seuratObj, 
                               libraryID = libraryID, 
                               studyID = studyID,
                               seuratFilter = FALSE)
  
  # Return some key pre-QC metrics 
  qcMetrics = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                "percent.ribo", "percent.hb")
  qcPlot = plotRNAFeatureList(seuratObj = seuratObj, 
                              features = qcMetrics)

  summaryDF = makeSumStatsDF(seuratObj = seuratObj)
  
  # Calculate library complexity scores and generate plot 
  complexityPlot = plotRNAcomplexity(seuratObj = seuratObj, cor = TRUE)
  
  # Define range of resolutions to test, create euclidean distance matrix and calculate silhouette coeffs
  message("--------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  clusteringRes = seq(minRes, maxRes, by = 0.1)
  prefiltSilcoeff = calcSilscores(seuratObj, clusteringRes)
  
  # Get mean sil_scores per resolution and select the maximum. Call our wrapper in utiliy
  maxSilRes = maxAvgSil(prefiltSilcoeff)
  message("The maximum avg sil score for prefiltered dataset is ", maxSilRes)
  
  # Now we use the max avg sil value to cluster data. Louvain clustering is probably fine
  # to identify doublets
  seuratObj =  FindClusters(seuratObj, resolution = maxSilRes)
  beforeQCclusters = DimPlot(seuratObj, reduction = "umap", label = TRUE)
  
  # Produce a list of whatever genes the user needs
  if (!is.null(queryGenes)) {
    beforeQCfeatures = plotGeneUMAP(seuratObj = seuratObj,
                                    queryGenes = queryGenes)
  }

  ######################################################################
  # STEP 3: Once we have defined our clusters, we can proceed to remove doublets with the wrapper
  # we wrote in the utility script. 
  message("------------------
  Finding doublets...
  -------------------")
  doubletIDs = scDblFinderClusters(seuratObj = seuratObj, 
                                     nrep = dblFindIter)
  message("There are ", length(doubletIDs), "from ", 
          dblFindIter, " iterations of scDblFinder")
  
  # First filtering step: remove consensus doublet IDs
  singlets = setdiff(colnames(countsMatrix), doubletIDs)
  message("There are ", length(singlets), 
          " cells remaining after removing doublets")
  
  # Add UMAP plot showing consensus doublets
  doubletsUMAP = DimPlot(seuratObj, cells.highlight = doubletIDs) + 
    custom_theme() + 
    ggtitle("Consensus doublets")
  
  # Add qc plots into list (Before filtering and decontamination)
  preQCmetrics = list(metricsPlot = qcPlot, 
                      metricsSummary = summaryDF,
                      libraryComplexity = complexityPlot, 
                      doubletsUMAP = doubletsUMAP)
  
  # Create a new seurat obj keeping only singlets
  seuratObj = CreateSeuratObject(counts=countsMatrix[, singlets],
                                  project = libraryID,
                                  min.cells = 10,
                                  min.features = 200)
  
  # Clean newly filtered counts matrix from ambient RNA using our own decontX wrapper
  # within the utility script. 
  message("---------------------------------------------------------
  Removing ambient RNA...because who doesn't like cleaner data? -
  ------------------------------------------------------------------")
  
  # This step will reduce the number of UMIs to less than the specified threshold above
  # because it's removing contaminating reads. This means we might have to harcode 
  # a lower bound threshold
  seuratObj = decontXremove(seuratObj = seuratObj)
  
  ##########################################################################
  # STEP4: Normalization and dimensionality reduction of cleaned counts matrix
  # In this case, the matrix will be further filtered according to nUMIs,
  # nFeatures and percentage of reads mapped to mitochondrial genome. 
  seuratObj = seuratSCTprocess(seuratObj = seuratObj,
                               libraryID = libraryID,
                               studyID = studyID,
                               minUMI = minUMI,
                               maxUMI = maxUMI,
                               minFeatures = minFeatures,
                               maxFeatures = maxFeatures,
                               maxMtPercent = maxMtPercent,
                               maxHbPercent = maxHbPercent,
                               arterialOrigin = arterialOrigin,
                               diseaseStatus = diseaseStatus,
                               sex = sex,
                               seuratFilter = TRUE,
                               setAutoThreshold = TRUE)
  
  # Add complexity scores into metadata
  nGenes = seuratObj$nFeature_RNA
  nUMIs = seuratObj$nCount_RNA
  seuratObj@meta.data$log10GenesPerUMI = log10(nGenes) / log10(nUMIs)
  
  # Before clustering we need to define the most optimal resolution... again
  message("-------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  
  filteredSilCoeff = calcSilscores(seuratObj = seuratObj,
                                   clusteringRes,
                                   clusteringAlg = ifelse(
                                     clusteringAlg == "louvain", 
                                     "louvain", 
                                     "leiden"))
  
  filteredMaxSil = maxAvgSil(filteredSilCoeff)
  
  message("The maximum avg sil score for the filtered data is ", filteredMaxSil)
  
  # Now we can cluster the data 
  message("Clustering data...")
  
  # Add arg to produce clusters using either Louvain or Leiden algorithm
  seuratObj = FindClusters(seuratObj,
                           resolution = filteredMaxSil,
                           algorithm = ifelse(clusteringAlg == "louvain", 1, 4))
  
  afterQCclusters = DimPlot(seuratObj, reduction = "umap", label = TRUE)
  
  # Plot post QC metrics for nUMIs, nGenes to compare with pre-filtered
  postQCmetrics = plotRNAFeatureList(seuratObj = seuratObj,
                                     features = qcMetrics)
  postQCSummary = makeSumStatsDF(seuratObj = seuratObj)
  postQCcomplexity = plotRNAcomplexity(seuratObj = seuratObj, cor = TRUE)
  postQCmetrics = list(metricsPlot = postQCmetrics,
                       metricsSummary = postQCSummary,
                       libraryComplexity = postQCcomplexity)
  
  # Produce list of gene expression plots after QC. 
  if (!is.null(queryGenes)) {
    afterQCfeatures = plotGeneUMAP(seuratObj = seuratObj,
                                   queryGenes = queryGenes)
  }
  
  message("---------------------------------------------------
  The input library has been succesfully processed :)
  ---------------------------------------------------")
  
  if (!is.null(queryGenes)) {
    Outputs = list(seuratObj = seuratObj,
                   doubletBarcodes = doubletIDs,
                   preQCplotsList = preQCmetrics,
                   postQCplotsList = postQCmetrics,
                   beforeQCclusters = beforeQCclusters,
                   afterQCclusters = afterQCclusters,
                   beforeQCfeatures = beforeQCfeatures,
                   afterQCfeatures = afterQCfeatures)
  } else {
    Outputs = list(seuratObj = seuratObj,
                   doubletBarcodes = doubletIDs,
                   preQCplotsList = preQCmetrics,
                   postQCplotsList = postQCmetrics,
                   beforeQCclusters = beforeQCclusters,
                   afterQCclusters = afterQCclusters)
  }
  
  return(Outputs)
  
}






