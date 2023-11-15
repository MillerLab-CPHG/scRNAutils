

#' Process scRNA libraries 
#' 
#' A function that will handle QC (doublet and ambient RNA removal) and processing for individual scRNA libraries
#' 
#' @param countsMatrix sparse count matrix required as initial input. As of now, make it accept 
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
  minRes = 0.2, 
  maxRes = 1.9, 
  dblFindIter = 3,
  setAutoThreshold = FALSE,
  regressMitoRibo = FALSE,
  rmMitoRiboVarGenes = FALSE,
  clusteringAlg = "louvain",
  queryGenes = NULL,
  minUMI = 500,
  maxUMI = 20000,
  minFeatures = 200,
  maxFeatures = 4000,
  maxMtPercent = 10,
  maxHbPercent = 1,
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
  preQCcols = as.numeric(dim(countsMatrix)[2])
  message(paste("Your prefiltered count matrix has", 
                dim(countsMatrix)[1], "genes and",
                preQCcols, "cells"))
  # QC (removal of doublets); Prelim round of clustering
  seuratObj = CreateSeuratObject(
    counts = countsMatrix,
    project = libraryID,
    min.cells = 10,
    min.features = 200
    )
  # Do a first round of clustering for doublet removal
  # Set the seurat_filter arg to FALSE
  seuratObj = seuratSCTprocess(
    seuratObj = seuratObj, 
    libraryID = libraryID, 
    studyID = studyID, 
    seuratFilter = FALSE
    )
  # Return some key pre-QC metrics 
  qcMetrics = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                "percent.ribo", "percent.hb")
  qcPlot = plotRNAFeatureList(
    seuratObj = seuratObj, 
    features = qcMetrics
    )
  summaryDF = makeSumStatsDF(seuratObj = seuratObj)
  complexityPlot = plotRNAcomplexity(seuratObj = seuratObj, cor = TRUE)
  
  # Define range of resolutions to test, create euclidean distance matrix and calculate silhouette coeffs
  message("--------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  clusteringRes = seq(minRes, maxRes, by = 0.1)
  prefiltSilcoeff = calcSilscores(seuratObj, clusteringRes)
  # Get mean sil_scores per resolution 
  if (!all(sapply(prefiltSilcoeff, is.null))) {
    clusterResPre = maxAvgSil(prefiltSilcoeff)
    message("The resolution with the maximum avg sil score for prefiltered dataset is ", 
            clusterResPre)
  } else {
    clusterResPre = minRes
  }
  # Use the max avg sil value to cluster data
  seuratObj =  FindClusters(seuratObj, resolution = clusterResPre)
  preQCclusters = DimPlot(seuratObj, reduction = "umap", label = TRUE)
  
  # Produce a list of whatever genes the user needs
  if (!is.null(queryGenes)) {
    preQCfeatures = plotGeneUMAP(
      seuratObj = seuratObj,
      queryGenes = queryGenes
      )
  }
  # STEP 3: Remove doublets with our scDblFinder wrapper
  message("------------------
  Finding doublets...
  -------------------")
  doubletIDs = scDblFinderClusters(
    seuratObj = seuratObj, 
    nrep = dblFindIter
    )
  message("There are ", length(doubletIDs), "from ", 
          dblFindIter, " iterations of scDblFinder")
  
  # First filtering step: remove consensus doublet IDs
  singlets = setdiff(colnames(countsMatrix), doubletIDs)
  message("There are ", length(singlets), 
          " cells remaining after removing doublets")
  doubletsUMAP = DimPlot(seuratObj, cells.highlight = doubletIDs) + 
    custom_theme() + 
    ggtitle("Consensus doublets")
  # Add pre QC plots 
  preQCmetrics = list(
    metricsPlot = qcPlot, 
    metricsSummary = summaryDF,
    libraryComplexity = complexityPlot, 
    doubletsUMAP = doubletsUMAP)
  # Create a new seurat obj keeping only singlets
  seuratObj = CreateSeuratObject(
    counts=countsMatrix[, singlets],
    project = libraryID,
    min.cells = 10,
    min.features = 200
    )
  message("---------------------------------------------------------
  Removing ambient RNA...because who doesn't like cleaner data? -
  ------------------------------------------------------------------")
  # This step will reduce the number of UMIs to less than the specified threshold above
  # because it's removing contaminating reads. This means we might have to harcode 
  # a lower bound threshold
  seuratObj = decontXremove(seuratObj = seuratObj)
  
  # STEP4: Normalization and dim reduction of cleaned counts matrix
  seuratObj = seuratSCTprocess(
    seuratObj = seuratObj,
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
    rmMitoRiboVarGenes = TRUE,
    seuratFilter = TRUE,
    setAutoThreshold = TRUE,
    regressMitoRibo = FALSE
    )
  # Add complexity scores into metadata
  nGenes = seuratObj$nFeature_RNA
  nUMIs = seuratObj$nCount_RNA
  seuratObj$log10GenesPerUMI = log10(nGenes) / log10(nUMIs)
  
  # Define the most optimal resolution... again
  message("-------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  filteredSilCoeff = calcSilscores(
    seuratObj = seuratObj,
    clusteringRes,
    clusteringAlg = ifelse(clusteringAlg == "louvain", 
                           "louvain", "leiden")
    )
  if (!all(sapply(prefiltSilcoeff, is.null))) {
    clusterResPos = maxAvgSil(filteredSilCoeff)
    message("The resolution with the maximum avg sil score for the filtered data is ", 
            clusterResPos)
  } else {
    clusterResPos = minRes
  }
  message("Clustering data...")
  # 1 = Louvain; 4 = Leiden
  seuratObj = FindClusters(
    seuratObj,
    resolution = clusterResPos,
    algorithm = ifelse(
      clusteringAlg == "louvain", 1, 4)
    )
  # Make QC plots
  message("Producing output plots...")
  postQCclusters = DimPlot(seuratObj, reduction = "umap", label = TRUE)
  postQCmetrics = plotRNAFeatureList(seuratObj = seuratObj,
                                     features = qcMetrics)
  postQCSummary = makeSumStatsDF(seuratObj = seuratObj)
  postQCcomplexity = plotRNAcomplexity(seuratObj = seuratObj, cor = TRUE)
  postQCcomplexityHist = plotRNAcomplexityHist(seuratObj = seuratObj)
  percentPlots = plotFeatureUMAPList(seuratObj = seuratObj,
    queryFeatures = c(
      "percent.ribo", 
      "percent.mt",
      "contamination_scores")
    )
  # Save post-QC plots and stats
  postQCcols = ncol(seuratObj)
  topGenes = plotTopAvgExpr(seuratObj = seuratObj)
  postQCmetrics = list(
    cellNumber = data.frame(
      preQCcells = preQCcols,
      postQCcells = postQCcols
    ),
    metricsPlot = postQCmetrics,
    metricsSummary = postQCSummary,
    libraryComplexity = postQCcomplexity,
    complexityHist = postQCcomplexityHist,
    contaminationScores = percentPlots$contamination_scores,
    percentRiboplot = percentPlots$percent.ribo,
    percentMTplot = percentPlots$percent.mt,
    topGenes = topGenes
    )
  # Produce list of gene expression plots after QC. 
  if (!is.null(queryGenes)) {
    postQCfeatures = plotGeneUMAP(
      seuratObj = seuratObj,
      queryGenes = queryGenes
      )
  }
  Outputs = list(
    seuratObj = seuratObj,
    doubletBC = doubletIDs,
    preQCstats = preQCmetrics,
    postQCstats = postQCmetrics,
    preQCclusters = preQCclusters,
    postQCclusters = postQCclusters)
  if (!is.null(queryGenes)) { 
    Outputs[["preQCfeatures"]] = preQCfeatures
    Outputs[["postQCfeatures"]] = postQCfeatures
  }
  message("---------------------------------------------------
  The input library has been succesfully processed :)
  ---------------------------------------------------")
  
  return(Outputs)
}






