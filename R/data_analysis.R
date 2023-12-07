


#' Identify and remove doublets. 
#' 
#' Internal helper function to automate doublet removal workflow and update metadata.
#' 
#' @param seuratObj: A Pre-processed seurat object with cell clusters annotations.
#' @param nrep A Numeric indicating the number of times to run the scDblFinder function.
#' @param multisampleDataset A Logical indicating whether the processed dataset has multiple samples
#' @param sample_ids if multisample_dataset = TRUE, A character that refers to the metadata column
#' that contains sample ids.
#' @return A vector of consensus doublet cell barcodes
scDblFinderClusters = function(
  seuratObj, 
  nrep, 
  multisampleDataset = FALSE, 
  sampleIDs = NULL
  ){
  # Create a temp single cell experiment object
  sampleSce = as.SingleCellExperiment(seuratObj)
  # If we have more than one sample in a dataset, use the samples parameter
  if (multisampleDataset) {
    sampleSceDbl = replicate(
      nrep, 
      scDblFinder::scDblFinder(
        sampleSce,
        samples = sampleIDs,
        clusters = "seurat_clusters"
        )
      )
  } else {
    sampleSceDbl = replicate(
      nrep, scDblFinder::scDblFinder(
        sampleSce, 
        clusters = "seurat_clusters"
        )
      )
  }
  newSeurat = lapply(sampleSceDbl, as.Seurat)
  # Get barcodes that were tagged as doublets
  res = lapply(newSeurat, 
               function(x){
                 x@meta.data %>% 
                   filter(scDblFinder.class == "doublet") %>% 
                   rownames()})
  rm(sampleSce)
  # Get consensus doublet calls between the desired number of runs
  consensusDbl = Reduce(intersect, res[1:length(res)])
  return(consensusDbl)
}


#' Removal of ambient RNA.
#' 
#' Internal decontX wrapper to standardize removal of ambient mRNA contamination
#' 
#' @param seuratObj: A seurat obj that has been filtered for doublets
#' @return A seurat object with a decontaminated raw counts matrix and 
#' a metadata column showing the contamination scores output by decontX. 
decontXremove = function(seuratObj) {
  
  # Run main decontX function
  res = celda::decontX(seuratObj@assays$RNA@counts)
  
  # Access decontaminated counts 
  # Decont matrix contains floating point numbers so might need to round them
  decontMatrix = res$decontXcounts
  
  # Add decont Matrix and contamination scores into Seurat obj
  seuratObj@assays$RNA@counts = decontMatrix
  seuratObj@meta.data$contamination_scores = res$contamination
  return(seuratObj)
}


#' Downstream processing of counts matrix
#' 
#' Internal helper function to standardize filtering, SCT normalization
#' and dimensionality reduction. This function will also add
#' some important metadata variables.  
#' 
#' @param seuratObj An input Seurat object.
#' @param seuratFilter A boolean indicating whether we wish to filter low quality cells. 
#' Set to FALSE in pre-processing round for doublet removal.
#' @param minUMI A numeric indicating min number of UMIs to keep.
#' @param maxUMI A numeric indicating max number of UMIs to keep.
#' @param minFeatures A numeric indicating the minimum number of genes cells should express to be included.
#' @param maxFeatures A numeric indicating the maximum number of genes cells should express to be included.
#' @param maxMtPercent A numeric indicating the max percentage of reads mapped to the mito genome. 
#' @param maxHbPercent A numeric indicating the max percentage of reads mapped to hemoglobin genes. 
#' @param libraryID A character indicating name of the sample.
#' @param studyID A character indicating name of the study. 
#' @param arterialOrigin A character indicating arterial bed of the library.
#' @param diseaseStatus A character indicating the disease status of the library (e.g., non-lesion, lesion).
#' @param sex A character string indicating sex of the subject. 
#' @param regressMito A boolean indicating whether to regress out Mitochondrial variance or not.
#' @param setAutoThreshold A boolean indicating whether we should do adative thresholding for removing
#' low quality cells. This approach uses the Median Absolute Deviation (MAD).  
#' @param rmMitoRiboVarGenes A boolean indicating whether we should remove Mito and Ribo genes from 
#' the set of highly variable features.  
#' @param autoSelectPCs A boolean indicating whether to determine the number of PCs in a data-driven way. 
#' 
#' @return A seurat object that has been SCT normalized, nearest neighbors graph and UMAP embeddings.
seuratSCTprocess = function(
  seuratObj, 
  libraryID,
  studyID,
  minUMI = 500,
  maxUMI = 20000,
  minFeatures = 200,
  maxFeatures = 4000,
  maxMtPercent = 10,
  maxHbPercent = 1,
  arterialOrigin = NULL, 
  diseaseStatus = NULL,
  sex = NULL,
  regressMitoRibo = FALSE,
  seuratFilter = FALSE,
  setAutoThreshold = TRUE,
  rmMitoRiboVarGenes = FALSE,
  autoSelectPCs = TRUE
  ){
  # Define sample and study IDs as core metadata values. 
  if (is.null(libraryID) || is.null(studyID)) { 
    msg = paste("Library and Study IDs are missing!")
    stop(msg)
    }
  # Add required metadata variables 
  seuratObj$sample = libraryID
  seuratObj$study = studyID
  
  # Vascular bed and sample disease status can be kept as optional metadata values.
  if (!is.null(arterialOrigin)) {
    seuratObj$arterialOrigin = arterialOrigin
  }
  if (!is.null(diseaseStatus)) {
    seuratObj$diseaseStatus = diseaseStatus
  }
  if (!is.null(sex)) { 
    seuratObj$sex = sex
  }
  # Quality control
  # Reads mapping to Mito genes
  seuratObj[["percent.mt"]] = PercentageFeatureSet(seuratObj, 
                                                    pattern = "^MT-")
  # Reads mapping to hemoglobin genes
  hbIndex = grep(rownames(seuratObj), pattern = "^HB[AB]")
  hbGenes = rownames(seuratObj)[hbIndex]
  seuratObj[["percent.hb"]] = PercentageFeatureSet(seuratObj, features = hbGenes)
  # Reads mapping to ribosomal genes
  riboIndex = grep(
    rownames(seuratObj), 
    pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"
    )
  riboGenes = rownames(seuratObj)[riboIndex]
  seuratObj[["percent.ribo"]] = PercentageFeatureSet(seuratObj,features = riboGenes)
  if (seuratFilter) {
    if (setAutoThreshold) {
      # Use the MAD to identify outliers based on joint QC metrics and use isOutlier()
      # default is 3 MADs below/above median, might need to try 2 MADs
      message("Setting MAD-based adaptive thresholds for cells filtering...")
      stats = cbind(
        log10(seuratObj$nCount_RNA),
        log10(seuratObj$nFeature_RNA),
        seuratObj$percent.mt,
        seuratObj$percent.ribo
        )
      outlying = adjOutlyingness(stats, only.outlyingness = TRUE)
      outliersIdx = isOutlier(outlying, type = "both", nmads = 1)
      outlierBarcodes = colnames(seuratObj)[outliersIdx]
      seuratObj = subset(
        seuratObj, 
        cells = outlierBarcodes, invert = TRUE
        )
      seuratObj = subset(
        seuratObj,
        subset = nFeature_RNA >= minFeatures & 
          nCount_RNA >= minUMI & 
          percent.hb <= maxHbPercent)
      seuratObj = subset(
        seuratObj,
        subset = contamination_scores <= 0.2
      )
    } else {
      seuratObj =  subset(
        seuratObj,
        subset =  nFeature_RNA <= maxFeatures & 
          nCount_RNA <= maxUMI & 
          percent.mt <= maxMtPercent & 
          percent.hb <= maxHbPercent
        )
      seuratObj = subset(
        seuratObj,
        subset = nFeature_RNA >= minFeatures 
        & nCount_RNA >= minUMI
        )
    }
  }
  # Calculate cell cycle scores
  Seurat::cc.genes.updated.2019
  s.genes = cc.genes.updated.2019$s.genes
  g2m.genes = cc.genes.updated.2019$g2m.genes
  seuratObj = CellCycleScoring(
    seuratObj, s.features = s.genes, 
    g2m.features = g2m.genes
    )
  # Normalize data, find variable genes, scale data and regress out cell cycle variance
  # SCT enables extraction of meaningful insights from more PCs so we'll set dims=1:30
  # By default, regress out cell cycle variance
  cellCycleVec = c("S.Score","G2M.Score")
  allCovariatesVec = c(cellCycleVec, "percent.mt", "percent.ribo")
  seuratObj = SCTransform(
    seuratObj,
    vst.flavor = "v2",
    vars.to.regress = if (regressMitoRibo) {
                            allCovariatesVec
                          } else {
                            cellCycleVec
                          }
  )
  # Remove mito and ribosomal genes from set of variable features
  if (rmMitoRiboVarGenes) {
    message("Removing Mito and Ribosomal genes from Variable features...")
    mitoIndex = grep(rownames(seuratObj), pattern = "^MT-")
    mitoGenes = rownames(seuratObj)[mitoIndex]
    varGenes = setdiff(
      VariableFeatures(seuratObj), 
      c(mitoGenes, riboGenes, "MALAT1")
      )
    VariableFeatures(seuratObj) = varGenes
  } 
  # Use VarGenes for dim reduction and clustering
  seuratObj = RunPCA(
    seuratObj, 
    features = VariableFeatures(seuratObj)
    ) 
  nPCs = selectPCs(seuratObj = seuratObj)
  message(paste0(nPCs, " PCs explain 90% of the variance in your library..."))
  seuratObj = FindNeighbors(
    seuratObj, 
    reduction = "pca", 
    dims = if(autoSelectPCs) 1:nPCs else 1:30, 
    k.param = 20
    ) %>%
    RunUMAP(
      dims = if (autoSelectPCs) 1:nPCs else 1:30, 
      n.neighbors = 30
      )
  return(seuratObj)
}

#' Helper function to make a summary df from metadata vars
#' 
#' This helper function will return summary stats for key QC metrics
#' as a dataframe.
#' 
#' @param seuratObj A Seurat object with nUMI, nFeatures and other QC metrics
#' @return A dataframe with summary statistics for several important QC metrics
#' 
#' @export
#' 
makeSumStatsDF = function(seuratObj) {
  metricsSummary = list(
    nUMIs = summary(seuratObj$nCount_RNA),
    nGenes = summary(seuratObj$nFeature_RNA),
    percentMT = summary(seuratObj$percent.mt),
    percentRibo = summary(seuratObj$percent.ribo),
    percentHb = summary(seuratObj$percent.hb))
  summaryDF = as.data.frame(do.call(rbind, metricsSummary))
  return(summaryDF)
}

#' Select n PCs
#'
#' This function will find the number of PCs explaining the largest
#' amount of variation in the data, which we need
#' for clustering and UMAP in a more data-driven way.
#' 
#' @param seuratObj A seurat obj with a computed DimReduc object
#' @return A numeric indicating the number of PCs we need to explain the largest variation
#'
selectPCs = function(seuratObj) { 
  # Extract eigen values and calc % variation
  message("Selecting the n PCs that capture 90% of variation...")
  eigenVals = seuratObj@reductions$pca@stdev
  explainedVar = eigenVals / sum(eigenVals) * 100
  cumVar = cumsum(explainedVar)
  co1 = which(cumVar > 90 & explainedVar < 5)[1]
  # Second metric: point where change in variation is less than 0.1%
  #varDelta = c(0, abs(diff(explainedVar)))
  #co2 = sort(which(varDelta > 0.1), decreasing = TRUE)[1]
  #nPCs = min(co1, co2)
  return(co1)
  }


#' Calculation of silhouette scores
#' 
#' Helper function to calculate clusters silhouette scores from a seurat object.
#' 
#' @param seuratObj A Seurat obj with reduced dims embedding like pca. Default now is set to PCA. 
#' @param clusteringRes A Numeric vector with clustering resolutions to test
#' @return A Data.frame with silhouette scores for a the provided clustering resolutions
#' @export
#' @examples 
#' \dontrun{
#' sil_scores = calc_sil_scores(seurat_obj, seq(0.3, 0.9, by=0.1))
#' }
calcSilscores = function(
  seuratObj, 
  clusteringRes,
  clusteringAlg = "louvain"
  ) {
  clusteringRes = clusteringRes
  silScoresList = vector("list", length(clusteringRes))
  message("Creating distance matrix from PCA embeddings...")
  distMatrix = stats::dist(
    x = Seurat::Embeddings(object = seuratObj, 
                           reduction = "pca")[, 1:30]
    )
  # Cluster data with each of the resolutions provided
  for (i in seq_along(clusteringRes)) { 
    tryCatch({
      message("Clustering data...")
      seuratObj = FindClusters(
        seuratObj,
        resolution = clusteringRes[i],
        algorithm = ifelse(clusteringAlg == "louvain", 1, 4)
        )
      clusters = seuratObj$seurat_clusters
      if (length(unique(clusters)) > 1) {
        message("Calculating silhouette scores for resolution: ", clusteringRes[i])
        sil = cluster::silhouette(
          x = as.numeric(x = as.factor(x = clusters)), 
          dist = distMatrix
          )
        silDf = data.frame(
          cluster=sil[, 1],
          neighbor=sil[, 2],
          silWidth=sil[, 3]
          )
        silScoresList[[i]] = silDf
        names(silScoresList)[[i]] = paste("res", as.character(clusteringRes[i]), sep = "_")
        }
      }, error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
      })
    }
  resolutionsSilDf = data.table::rbindlist(
    silScoresList,
    idcol = TRUE)
  return(resolutionsSilDf)
}


#' Internal helper function to get the max sil score from calc_sil_scores()
#' 
#' @param silScoresDf A data.frame output by calc_sil_scores()
#' @return A numeric indicating the maximum average sil sclore value across tested resolutions
maxAvgSil = function(silScoresDf) 
  {
  avgSilScores = silScoresDf %>%
    group_by(.id) %>%
    summarize(meanSil = mean(silWidth))
  maxSilRes = as.character(avgSilScores[which.max(avgSilScores$meanSil), 1])
  maxSilRes = as.numeric(str_split(maxSilRes, "_")[[1]][2])
  return(maxSilRes)
}


#' Correlation of expression for a target gene across the entire transcriptome. 
#'
#' This function will calculate the correlation between a target gene and other genes in a specific group of cells. 
#' 
#' @param seuratObj A Seurat object.
#' @param assay A character vector indicating the assay from which to get the counts.
#' @param targetGene A string indicating the target gene to correlate with the rest of the transcriptome. 
#' @param corMethod A string indicating the correlation method to use, either spearman or pearson
#' @param annoVariable A character vector indicating the name of the metadata column that has the cell type annotations
#' @param cellTypes A character vector indicating the cell type annotations in which to calculate the correlations
#' 
#' @return A data.frame with the correlation of the target gene with all other genes in the count matrix. 
#' @export
#' 
calcGeneCors = function(
  seuratObj,
  assay = "SCT",
  targetGene, 
  corMethod = "pearson",
  annoVariable = NULL,
  cellTypes = NULL
  ) {
  # Set default assay
  DefaultAssay(seuratObj) = assay
  # Create a vector with the names of the cells of interest
  metadataDf = seuratObj@meta.data
  if (!is.null(annoVariable)) {
    if (is.null(cellTypes)) {
      stop("Please provide a list of cell type annotations in ", 
           annoVariable)
    }
    cellsVector = metadataDf %>%
      dplyr::filter(annoVariable %in% cellTypes)
  } else {
    cellsVector = colnames(seuratObj)
  }
  # Process matrix to run correlations
  if (!is.null(annoVariable) && !is.null(cellTypes)) {
    geneExpression = seuratObj@assays$SCT@data[, cellsVector]
  } else {
    geneExpression = seuratObj@assays$SCT@data
  }
  reshapedMatrix = t(as.matrix(geneExpression))
  targetGeneVec = as.numeric(reshapedMatrix[, targetGene])
  corsList = apply(
    reshapedMatrix, 2, 
    function(x){
      cor.test(
        targetGeneVec, 
        x,
        method = corMethod
        )
      }
    )
  corsListTidy = lapply(corsList, broom::tidy)
  corsDf =  data.table::rbindlist(corsListTidy)
  geneNames = names(corsList)
  corsDfGenes = corsDf %>%
    mutate(gene=geneNames,
           FDR = p.adjust(p.value, "BH")) %>%
    arrange(desc(estimate)) %>%
    drop_na() %>%
    mutate(geneIndex = seq_along(gene)) %>%
    filter(!gene == targetGene) %>%
    dplyr::rename(corEstimate = estimate)
  return(corsDfGenes)
}










