


#' Tidy Seurat counts into a data.frame with cell types of interest 
#' 
#' @param seurat_obj Seurat object from which we'll extract the normalized counts.
#' @param assay A character indicating the name of the assay from which to extract the counts.
#' This can be set to SCT or RNA. 
#' @param target_coldata A character indicating the name of the column with the cluster IDs or cell type annotations.
#' @param target_anno A character vector with the cluster IDs or names of cell types that will be included in the df. 
#' 
#' @return A data.frame of cells x genes with corresponding normalized counts. This df also includes a column with the 
#' cell barcode matching the cluster ID or cell type annotation.  
#' @export
#' 
tidySeuratCounts = function(
  seuratObj, 
  assay="SCT", 
  targetColData=NULL, 
  targetAnno=NULL
  ) {
  # Get counts matrix
  if (assay == "SCT") {
    normCounts = seuratObj[["SCT"]]$data
  } else if (assay == "RNA") {
    normCounts = seuratObj[["RNA"]]$data
  } else {
    msg = ("The provided assay should be either SCT or RNA!")
    stop(msg)
  }
  normMatrixDf = as.data.frame(t(as.matrix(normCounts)))
  if (!is.null(targetAnno) && !is.null(targetColData)) {
    metaDf = seuratObj[[]]
    # Subset df and name barcodes with cell type annotations
    metaDfsubset = metaDf[metaDf[[targetColData]] %in% targetAnno, ]
    targetBarcodes = rownames(metaDfsubset)
    names(targetBarcodes) = metaDfsubset[[targetColData]]
    
    # Subset counts df to keep only cell types of interest
    countsDfsubset = normMatrixDf[targetBarcodes,]
    countsDfsubset$annotation = names(targetBarcodes)
    return(countsDfsubset)
  } else {
    return(normMatrixDf)
  }
} 

