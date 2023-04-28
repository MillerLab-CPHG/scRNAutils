


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
tidy_seurat_counts = function(seurat_obj, assay="SCT", 
                              target_coldata=NULL, 
                              target_anno=NULL) {
  
  # Get counts matrix
  if (assay=="SCT") {
    norm_counts = seurat_obj@assays$SCT@data
  } else if (assay=="RNA") {
    norm_counts = seurat_obj@assays$RNA@data
  } else {
    msg = ("The provided assay should be either SCT or RNA!")
    stop(msg)
  }
  
  norm_matrix_counts_df = as.data.frame(t(as.matrix(norm_counts)))
  
  # What if we want to extract the entire counts matrix
  if (!is.null(target_anno) & !is.null(target_coldata)) {
    meta_df = seurat_obj@meta.data
    
    # Subset df
    meta_df_subset = meta_df[meta_df[[target_coldata]] %in% target_anno, ]
    target_barcodes = rownames(meta_df_subset)
    names(target_barcodes) = meta_df_subset[[target_coldata]]
    
    # Subset counts df to keep only cell types of interest
    counts_df_subset = norm_matrix_counts_df[target_barcodes,]
    
    # Add a column matching the cell barcode to its corresponding annotation 
    counts_df_subset$annotation = names(target_barcodes)
    return(counts_df_subset)
  } else {
    return(norm_matrix_counts_df)
  }
} 

