library(Seurat)
library(tidyverse)
library(data.table)
library(scDblFinder)
library(celda)
library(sctransform)
library(cluster)


# For now we need to source the utily script to load our own helper functions
# source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/project_scripts/scRNAutils/R/utility.R")
# source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/project_scripts/scRNAutils/R/data_analysis.R")
# source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/project_scripts/scRNAutils/R/data_manipulation.R")
# source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/project_scripts/scRNAutils/R/visualization.R")


#' Process scRNA libraries 
#' 
#' A function that will handle QC (doublet and ambient RNA removal) and processing for individual scRNA libraries
#' 
#' @param count_matrix sparse count matrix required as initial input. As of now, make it accept 
#' a filtered_feature_bc_matrix with a Seurat::Read10X or a txt file. Maybe this is better
#' handled out of the function since we'd also be dealing with h5 files. Maybe we can do this later
#' but not priority right now
#' @param library_id A character indicating the id for the library
#' @param min_res A numeric indicating lower boundary for sil analysis
#' @param max_res A numeric indicating upper boundary for sil analysis
#' @param dbl_remove_iter A numeric indicating how many iterations we need to run for getting consensus doublets
#' @param genes_of_interest A character vector with genes of interest for checking expression before and after QC. 
#' @return A list where the first object is the processed Seurat obj and remaining objects are UMAP plots before and after QC. 
#' @export
#' @examples 
#' seurat_outs = doitall(counts_matrix, "rpe004", "pan_et_al")

doitall = function(counts_matrix, library_id, study_id,
                   min_res=0.3, max_res=1.9, dbl_remove_iter=3,
                   genes_of_interest=NULL,
                   min_UMI=200 ,
                   max_UMI=20000,
                   min_features=200,
                   max_features=4000,
                   max_mt_percent=10,
                   max_hb_percent=5) {
  
  # Validate counts matrix input
  if (!is(counts_matrix, "dgCMatrix")) {
    if(is(counts_matrix, "Matrix")) {
      counts_matrix = as.sparse(counts_matrix)
    } else {
      stop("Your counts matrix should be provided in matrix or sparse matrix format (dgCMatrix)!. Got ", 
           class(counts_matrix))
    }
  }
  
  # How many rows and columns do this matrix have?
  message(paste("Your prefiltered count matrix has", dim(counts_matrix)[1], "rows and",
                dim(counts_matrix)[2], "columns"))
  
  # Load genes for regressing out cell cycle variance
  cc.genes.updated.2019
  s.genes = cc.genes.updated.2019$s.genes
  g2m.genes = cc.genes.updated.2019$g2m.genes
  
  # QC (removal of doublets)
  # To do this we need to create a Seurat object and do a first round of clustering
  # Keep genes expressed in at least 10 cells and cells with more than 200 expressed genes
  seurat_obj = Seurat::CreateSeuratObject(counts = counts_matrix,
                                          project = library_id,
                                          min.cells = 10,
                                          min.features = 200)
  
  # Do a first round of clustering for doublet removal
  # Set the seurat_filter arg to FALSE. We dont need to remove low quality cells
  # for this first round of clustering 
  seurat_obj = Seurat_SCT_process(seurat_obj = seurat_obj, 
                                  seurat_filter = FALSE,
                                  sample_id = library_id, 
                                  study_name = study_id,
                                  min_UMI = min_UMI,
                                  max_UMI = max_UMI,
                                  min_features = min_features,
                                  max_features = max_features,
                                  max_mt_percent = max_mt_percent,
                                  max_hb_percent = max_hb_percent)
  
  # How do we define the right clustering parameters? Maybe we need to calc sil scores
  # Define range of resolutions to test, create euclidean distance matrix and calculate silhouette coeffs
  # Min and max resolutions default to 0.3 and 1.9 but can be overwritten to any values really
  message("--------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  clustering_res = seq(min_res, max_res, by=0.1)
  prefiltered_sil_scores = calc_sil_scores(seurat_obj, clustering_res)
  
  # Get mean sil_scores per resolution and select the maximum. Call our wrapper in utiliy
  max_sil_res = max_avg_sil_score(prefiltered_sil_scores)
  message("The maximum avg sil score for prefiltered dataset is ", max_sil_res)
  
  
  # Now we use the max avg sil value to cluster data 
  seurat_obj =  FindClusters(seurat_obj, resolution = max_sil_res)
  before_qc_clusters = DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  
  # Produce a list of whatever genes the user needs
  if (!is.null(genes_of_interest)) {
    before_qc_features = seurat_gene_plot_list(seurat_object = seurat_obj,
                                               genes_of_interest = genes_of_interest)
  }

  ######################################################################
  # STEP 3: Once we have defined our clusters, we can proceed to remove doublets with the wrapper
  # we wrote in the utility script. Consensus doublets are obtained from however many iterations
  # the user defines. Default number of iterations is 3 but can set to whatever one desires. 
  message("------------------
  Finding doublets...
  -------------------")
  doublet_ids = scDblFinder_clusters(seurat_obj = seurat_obj, 
                                     nrep = dbl_remove_iter )
  message("There are ", length(doublet_ids), "from ", 
          dbl_remove_iter, " iterations of scDblFinder")
  
  # First filtering step: remove consensus doublet IDs
  singlets = setdiff(colnames(counts_matrix), doublet_ids)
  message("There are ", length(singlets), " cells remaining after removing doublets")
  
  # Create a new seurat obj keeping only singlets
  seurat_obj = CreateSeuratObject(counts = counts_matrix[, singlets],
                                  project = library_id,
                                  min.cells = 10,
                                  min.features = 200)
  
  # Clean newly filtered counts matrix from ambient RNA using our own decontX wrapper
  # within the utility script. Our wrapper returns a seurat object with a decontaminated
  # raw counts matrix
  message("---------------------------------------------------------
  Removing ambient RNA...because who doesn't like cleaner data? -
  ------------------------------------------------------------------")
  
  seurat_obj = decontX_remove(seurat_obj = seurat_obj)
  
  ##########################################################################
  # STEP4: Normalization and dimensionality reduction of cleaned counts matrix
  # In this case, the matrix will be further filtered accroding to nUMIs,
  # nFeatures and percentage of reads mapped to mitochondrial genome. 
  seurat_obj = Seurat_SCT_process(seurat_obj = seurat_obj,
                                  seurat_filter = TRUE,
                                  sample_id = library_id,
                                  study_name = study_id,
                                  min_UMI = min_UMI,
                                  max_UMI = max_UMI,
                                  min_features = min_features,
                                  max_features = max_features,
                                  max_mt_percent = max_mt_percent,
                                  max_hb_percent = max_hb_percent)
  
  # Before clustering we need to define the most optimal resolution... again
  message("-------------------------------------------------------------------------------
  Calculating maximum avg silhouette score for the provided clustering resolutions
  --------------------------------------------------------------------------------")
  filtered_sil_scores = calc_sil_scores(seurat_obj, clustering_res)
  filtered_max_sil_score = max_avg_sil_score(filtered_sil_scores)
  
  message("The maximum avg sil score for the filtered data is ", filtered_max_sil_score)
  
  # Now we can cluster the data 
  message("Clustering data...")
  seurat_obj = FindClusters(seurat_obj, resolution = filtered_max_sil_score)
  after_qc_clusters = DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  
  # Produce list of gene expression plots after QC. 
  if (!is.null(genes_of_interest)) {
    after_qc_features = seurat_gene_plot_list(seurat_object = seurat_obj,
                                              genes_of_interest = genes_of_interest)
  }
  
  message("---------------------------------------------------
  The input library has been succesfully processed :)
  ---------------------------------------------------")
  seurat_outputs = list(seurat_obj=seurat_obj, 
                        doublet_barcodes = doublet_ids,
                        before_qc_clusters=before_qc_clusters,
                        after_qc_clusters=after_qc_clusters,
                        before_qc_features = before_qc_features,
                        after_qc_features = after_qc_features)
  
  return(seurat_outputs)
  
}






