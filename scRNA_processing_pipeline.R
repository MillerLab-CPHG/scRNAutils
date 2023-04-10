library(Seurat)
library(tidyverse)
library(data.table)
library(sctransform)
library(cluster)


# Write a function that will do all of the scRNA-seq QC and processing
# for each sequencing library. We should make it as flexible as possible

# What are the inputs we should use? Start with a count matrix or a h5 file. 
# Args count_matrix count matrix required as initial input. As of now, make it accept 
# a filtered_feature_bc_matrix with a Seurat::Read10X or a txt file. Maybe this is better
# handled out of the function since we'd also be dealing with h5 files. Maybe we can do this later
# but not priority right now
# Args library_id A character indicating the id for the library


doitall = function(counts_matrix, library_id, study_id) {
  
  # Load data 
  if (!is(counts_matrix, "dgCMatrix")) {
    stop("Your counts matrix should be provided in sparse format (dgCMatrix)!")
  }
  
  # How many rows and columns do this matrix have?
  message(paste("Your prefiltered count matrix has", dim(counts_matrix)[1], "rows and",
                dim(counts_matrix)[2], "columns"))
  
  # QC (removal of doublets)
  # To do this we need to create a Seurat object and do a first round of clustering
  # Keep genes expressed in at least 10 cells and cells with more than 200 expressed genes
  seurat_obj = Seurat::CreateSeuratObject(counts = counts_matrix,
                                          project = library_id,
                                          min.cells = 10
                                          min.features = 200)
  
  # Do a first round of clustering for doublet removal
  # Set the seurat_filter arg to FALSE. We dont need to remove low quality cells
  # for this first round of clustering 
  seurat_obj = Seurat_SCT_process(rpe004_seurat_sct, 
                                  seurat_filter = FALSE,
                                  sample_id = library_id, 
                                  study_name = study_id)
  
  # How do we define the right clustering parameters? Maybe we need to calc sil scores
  # to determine best parameters before and after removing doublets. This might be important
  # since the actual resolution will depend on the number of cells in the library
  
  # Define range of resolutions to test and get PCA embeddings
  clustering_res = seq(0.3, 1.9, by=0.1)
  prefiltered_sil_scores = calc_sil_scores(seurat_obj, clustering_res)
  
  # Get mean sil_scores per resolution and select the maximum
  # Need to define the actual colname for the resolution vars
  avg_sil_scores = prefiltered_sil_scores %>%
    group_by(res) %>%
    summarize(mean_sil = mean(sil_width))
  max_sil_res = avg_sil_scores[which.max(avg_sil_scores$mean_sil), 1]
  

}






