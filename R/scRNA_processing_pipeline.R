

#' Process scRNA libraries. 
#' 
#' A function that will handle QC (doublet and ambient RNA removal) and robust processing for individual scRNA libraries.
#' 
#' @param countsMatrix sparse count matrix required as initial input. As of now, make it accept 
#'  a filtered_feature_bc_matrix with a Seurat::Read10X or a txt file. Maybe this is better
#'  handled out of the function since we'd also be dealing with h5 files. Maybe we can do this later
#'  but not priority right now
#' @param minRes A numeric indicating lower boundary for silhouette analysis.
#' @param maxRes A numeric indicating upper boundary for silhouette analysis.
#' @param dblFindIter A numeric indicating how many iterations we need to run for getting consensus doublets.
#' @param clusteringAlg A character vector indicating whether to perform louvain or leiden clustering.
#' @param queryFeatures A character vector with genes of interest for checking expression before and after QC.
#' @param UMI A boolean indicating whether the read type is a UMI.
#' @param makeAnnData A boolean indicating whether we want to generate an .h5ad file for the processed Seurat object.
#' @param annDataParentDir A character vector indicating the parent directory where we want to save the produced .h5ad object.
#' @inheritParams seuratSCTprocess
#' @return A list where the first object is the processed Seurat obj and a bunch of other QC stats and plots. 
#' @export
#' 
doItAll = function(
  countsMatrix, 
  libraryID,
  studyID,
  tissueSource,
  seqPlatform = "10x",
  UMI = TRUE,
  dblFindIter = 3,
  seuratFilter = TRUE,
  minRes = 0.2,
  maxRes = 1.9,
  clusteringAlg = "louvain",
  makeAnnData = FALSE,
  annDataParentDir = NULL,
  queryFeatures = NULL,
  ...
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
  
  # Do a first round of clustering for doublet removal.
  seuratObj = CreateSeuratObject(
    counts = countsMatrix,
    project = libraryID,
    min.cells = 10,
    min.features = 200
    )
  seuratObj = seuratSCTprocess(
    seuratObj = seuratObj, 
    libraryID = libraryID, 
    studyID = studyID, 
    tissueSource = tissueSource,
    seuratFilter = FALSE,
    regressMitoRibo = FALSE,
    rmMitoRiboVarGenes = FALSE
    )
  # Return some key pre-QC metrics. 
  qcMetrics = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                "percent.ribo", "percent.hb")
  qcPlot = plotRNAFeatureList(
    seuratObj = seuratObj, 
    features = qcMetrics
    )
  summaryDF = makeSumStatsDF(seuratObj = seuratObj, seqPlatform = seqPlatform)
  complexityPlot = plotRNAcomplexity(seuratObj = seuratObj, cor = TRUE, UMI = UMI)
  
  # Define res range, create euclidean dist matrix and calculate silhouette coeffs.
  message("--------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  clusteringRes = seq(minRes, maxRes, by = 0.1)
  prefiltSilcoeff = calcSilScores(seuratObj, clusteringRes)
  
  # Get mean sil_scores per resolution. 
  if (!all(sapply(prefiltSilcoeff, is.null))) {
    clusterResPre = maxAvgSil(prefiltSilcoeff)
    message("The resolution with the maximum avg sil score for prefiltered dataset is ", 
            clusterResPre)
  } else {
    clusterResPre = minRes
  }
  
  # Use the max avg sil value to cluster data.
  seuratObj =  FindClusters(seuratObj, resolution = clusterResPre)
  preQCclusters = DimPlot(seuratObj, reduction = "umap", label = TRUE)
  
  # Produce a list of whatever genes the user needs.
  if (!is.null(queryFeatures)) {
    preQCfeatures = plotFeatureUMAPList(
      seuratObj = seuratObj,
      queryFeatures =  queryFeatures
      )
  }
  
  # STEP 3: Remove doublets with our scDblFinder wrapper.
  message("------------------
  Finding doublets...
  -------------------")
  doubletIDs = scDblFinderClusters(
    seuratObj = seuratObj, 
    nrep = dblFindIter
    )
  message("There are ", length(doubletIDs), "from ", 
          dblFindIter, " iterations of scDblFinder")
  
  # First filtering step: remove consensus doublet IDs.
  singlets = setdiff(colnames(countsMatrix), doubletIDs)
  message("There are ", length(singlets), 
          " cells remaining after removing doublets")
  doubletsUMAP = DimPlot(seuratObj, cells.highlight = doubletIDs) + 
    custom_theme() + 
    ggtitle("Consensus doublets")
  
  # Create a new seurat obj keeping only singlets.
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
  
  contScores = plotRNAFeature(
    seuratObj = seuratObj, 
    feature = "contamination_scores"
    )
  
  # Add pre QC plots. 
  preQCmetrics = list(
    metricsPlot = qcPlot, 
    metricsSummary = summaryDF,
    libraryComplexity = complexityPlot, 
    doubletsUMAP = doubletsUMAP,
    contaminationScores = contScores
  )
  
  # STEP 4: Normalization and dim reduction of cleaned counts matrix.
  # Args inherited from SCTprocess.
  seuratObj = seuratSCTprocess(
    seuratObj = seuratObj,
    libraryID = libraryID,
    studyID = studyID,
    tissueSource = tissueSource,
    seuratFilter = seuratFilter,
    seqPlatform = seqPlatform,
    ...
    )
  
  # Add complexity scores into metadata.
  nGenes = seuratObj$nFeature_RNA
  nReads = seuratObj$nCount_RNA
  seuratObj$log10GenesPerReads = log10(nGenes) / log10(nReads)

  # Define the most optimal resolution.
  message("-------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  filteredSilCoeff = calcSilScores(
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
    algorithm = ifelse(clusteringAlg == "louvain", 1, 4))
  
  # Make QC plots.
  message("Producing output plots...")
  postQCclusters = DimPlot(seuratObj, reduction = "umap", label = TRUE)
  percentPlots = plotFeatureUMAPList(seuratObj = seuratObj,
    queryFeatures = c(
      "percent.ribo", 
      "percent.mt",
      "contamination_scores")
    )
  
  # Save post-QC plots and stats.
  postQCcols = ncol(seuratObj)
  postQCmetrics = list(
    cellNumber = data.frame(preQCcells = preQCcols, postQCcells = postQCcols),
    metricsPlot = plotRNAFeatureList(
      seuratObj = seuratObj, 
      features = c(qcMetrics, "contamination_scores")
      ),
    metricsSummary = makeSumStatsDF(seuratObj = seuratObj, seqPlatform = seqPlatform),
    libraryComplexity = plotRNAcomplexity(
      seuratObj = seuratObj, 
      cor = TRUE , 
      UMI = UMI
      ),
    complexityHist = plotRNAcomplexityHist(seuratObj = seuratObj),
    contaminationScores = percentPlots$contamination_scores,
    percentRiboplot = percentPlots$percent.ribo,
    percentMTplot = percentPlots$percent.mt,
    topGenes = plotTopAvgExpr(seuratObj = seuratObj, makeBoxplot = TRUE)
    )
  
  # Produce list of gene expression plots after QC. 
  if (!is.null(queryFeatures)) {
    postQCfeatures = plotFeatureUMAPList(
      seuratObj = seuratObj,
      queryFeatures = queryFeatures
      )
  }
  # Convert to Anndata.
  if (makeAnnData) {
    if (!dir.exists(annDataParentDir)) {
      message("Directory for output .h5ad file doesn't exist, creating one...")
      dir.create(annDataParentDir)
    }
    message("Converting Seurat object into AnnData object...")
    tryCatch({
      # Need raw counts for scvi-tools modeling.
      assays = c("SCT", "RNA")
      for (assay in assays) {
        if (assay == "RNA") {
          #Convert RNA assay to Assay3.
          seuratObj[["RNA"]] = suppressWarnings(
            as(object = seuratObj[["RNA"]], Class = "Assay")
          )
        }
        sceasy::convertFormat(
          obj = seuratObj,
          from = c("seurat"),
          to = c("anndata"),
          assay = assay,
          main_layer = "counts",
          drop_single_values = FALSE,
          outFile = paste0(
            annDataParentDir, "/", 
            unique(seuratObj$sample), "_", assay,
            "_counts.h5ad"
            )
          )
        if (assay == "RNA") {
          #Convert RNA assay back to Assay5.
          seuratObj[["RNA"]] = suppressWarnings(
            as(object = seuratObj[["RNA"]], Class = "Assay5")
          )
        }
      }
    }, 
      error = function(e) {
        cat("Sorry the object couldn't be converted:", 
            conditionMessage(e))
        }
    )
  }
  
  outputs = list(
    seuratObj = seuratObj,
    doubletBC = doubletIDs,
    preQCstats = preQCmetrics,
    postQCstats = postQCmetrics,
    preQCclusters = preQCclusters,
    postQCclusters = postQCclusters
    )
  if (!is.null(queryFeatures)) { 
    outputs[["preQCfeatures"]] = preQCfeatures
    outputs[["postQCfeatures"]] = postQCfeatures
  }
  
  message("---------------------------------------------------
  The input library has been succesfully processed :)
  ---------------------------------------------------")
  
  return(outputs)
}






