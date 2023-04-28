


#' Identify and remove doublets. 
#' 
#' Internal helper function to automate doublet removal workflow and update metadata.
#' 
#' @param seurat_obj: A Pre-processed seurat object with cell clusters annotations.
#' @param nrep A Numeric indicating the number of times to run the scDblFinder function.
#' @param multisample_dataset A Logical indicating whether the processed dataset has multiple samples
#' @param sample_ids if multisample_dataset = TRUE, A character that refers to the metadata column
#' that contains sample ids.
#' @return A vector of consensus doublet cell barcodes
scDblFinder_clusters = function(seurat_obj, nrep, 
                                multisample_dataset=FALSE, 
                                sample_ids=NULL) {
  
  sample_sce = as.SingleCellExperiment(seurat_obj)
  
  # If we have more than one sample in a dataset, use the samples parameter
  if (multisample_dataset) {
    sample_sce_dbl = replicate(nrep, scDblFinder::scDblFinder(sample_sce,
                                                 samples = sample_ids,
                                                 clusters = "seurat_clusters"))
  } else {
    sample_sce_dbl = replicate(nrep, scDblFinder::scDblFinder(sample_sce, 
                                                 clusters = "seurat_clusters"))
  }
  # Convert back to seurat objs
  new_seurat = lapply(sample_sce_dbl, as.Seurat)
  res = lapply(new_seurat, 
               function(x){x@meta.data %>% filter(scDblFinder.class=="doublet") %>% rownames()})
  
  # Get consensus doublet calls between the desired number of runs
  consensus_dbl_ids = Reduce(intersect, res[1:length(res)])
  return(consensus_dbl_ids)
}


#' Removal of ambient RNA
#' 
#' Internal decontX wrapper to standardize removal of ambient mRNA contamination
#' 
#' @param seurat_obj: A seurat obj that has been filtered for doublets
#' @return A seurat object with a decontaminated raw counts matrix
decontX_remove = function(seurat_obj) {
  decontX_results = celda::decontX(seurat_obj@assays$RNA@counts)
  decontaminated_matrix = decontX_results$decontXcounts
  seurat_obj@assays$RNA@counts = decontaminated_matrix
  return(seurat_obj)
}


#' Downstream processing of counts matrix
#' 
#' Internal helper function to standardize filtering, SCT normalization
#' and dimensionality reduction. This function will also add
#' some important metadata variables.  
#' 
#' @param seurat_obj An input seurat object.
#' @param seurat_filter A boolean indicating whether do we wish to filter. 
#' Set to FALSE in pre-processing round for doublet removal.
#' @param min_UMI A numeric indicating min number of UMIs to keep.
#' @param max_UMI A numeric indicating max number of UMIs to keep.
#' @param min_features A numeric indicating the minimum number of genes cells should express to be included.
#' @param max_features A numeric indicating the maximum number of genes cells should express to be included.
#' @param max_mt_percent A numeric indicating the max percentage of reads mapped to the mito genome. 
#' @param max_hb_percent A numeric indicating the max percentage of reads mapped to hemoglobin genes. 
#' @param library_id A character indicating name of the sample.
#' @param study_name A character indicating name of the study. 
#' @param arterial_origin A character indicating arterial bed of the library.
#' @param disease_status A character indicating the disease status of the library (e.g., non-lesion, lesion).
#' @return A seurat object that has been SCT normalized, nearest neighbors graph and UMAP embeddings.
Seurat_SCT_process = function(seurat_obj, 
                              seurat_filter=FALSE,
                              min_UMI=200 ,
                              max_UMI=20000,
                              min_features=200,
                              max_features=4000,
                              max_mt_percent=10,
                              max_hb_percent=5,
                              library_id, 
                              study_name, 
                              arterial_origin=NULL, 
                              sample_disease_status=NULL,
                              sex=NULL){
  
  # Define sample and study IDs as core metadata values. 
  seurat_obj$sample = library_id
  seurat_obj$study = study_name
  
  # Vascular bed and sample disease status can be kept as optional metadata values.
  if (!is.null(arterial_origin)) {
    seurat_obj$arterial_origin = arterial_origin
  }
  
  if (!is.null(sample_disease_status)) {
    seurat_obj$sample_disease_status = sample_disease_status
  }
  
  if (!is.null(sex)) { 
    seurat_obj$sex = sex
  }
  
  # Check for mt, hb percentage and other quality metrics 
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, 
                                                    pattern = "^MT-")
  
  # Define hemoglobin genes
  hb_index = grep(rownames(seurat_obj), pattern = "^HB[AB]")
  hb_genes = rownames(seurat_obj)[hb_index]
  seurat_obj[["percent.hb"]] = PercentageFeatureSet(seurat_obj, 
                                                    features = hb_genes)
  
  # Filter low quality cells (start with parameters used in the paper)
  # NOTE: Skip this step during the first processing round to be used for doublet detection
  if (seurat_filter) {
    seurat_obj =  subset(seurat_obj, 
                         subset = nFeature_RNA >= min_features & nFeature_RNA <= max_features & nCount_RNA >= min_UMI & nCount_RNA <= max_UMI & percent.mt <= max_mt_percent & percent.hb <= max_hb_percent)
  }
  
  # Calculate cell cycle scores
  seurat_obj = CellCycleScoring(seurat_obj, s.features = s.genes, 
                                g2m.features = g2m.genes)
  
  # Normalize data, find variable genes, scale data and regress out cell cycle variance
  # SCT enables extraction of meaningful insights from more PCs so we'll set dims=1:30
  # SCTransform Arg vst.flavor="v2" internally uses glmGamPoi
  seurat_obj = SCTransform(seurat_obj, vst.flavor = "v2", 
                           vars.to.regress = c("S.Score","G2M.Score")) %>%
    RunPCA() %>% 
    FindNeighbors(reduction = "pca", dims = 1:30, k.param = 20) %>%
    RunUMAP(dims = 1:30, n.neighbors = 30)
  
  return(seurat_obj)
}


#' Calculation of sihouette scores
#' 
#' Helper function to calculate clusters silhouette scores from a seurat object.
#' 
#' @param seurat_obj A Seurat obj with reduced dims embedding like pca. Default now is set to PCA. 
#' @param clustering_res A Numeric vector with clustering resolutions to test
#' @return A Data.frame with silhouette scores for a the provided clustering resolutions
#' @export
#' @examples 
#' \dontrun{
#' sil_scores = calc_sil_scores(seurat_obj, seq(0.3, 0.9, by=0.1))
#' }
calc_sil_scores = function(seurat_obj, 
                           clustering_res) {
  clustering_res = clustering_res
  sil_scores_list = list()
  
  # Create distance matrix from seurat obj reduced dims embedding 
  message("Creating distance matrix from PCA embeddings...")
  dist_matrix = stats::dist(x = Seurat::Embeddings(object = seurat_obj[["pca"]])[, 1:30])
  
  # Cluster data with each of the resolutions provided
  for (i in seq_along(clustering_res)) { 
    seurat_obj = FindClusters(seurat_obj, 
                              resolution = clustering_res[i])
    clusters = seurat_obj$seurat_clusters
    message("Calculating silhouette scores for defined resolutions...")
    
    sil = cluster::silhouette(x = as.numeric(x = as.factor(x = clusters)), 
                     dist = dist_matrix)
    # Add silhouette scores back into seurat obj metadata
    #seurat_obj$silhouette_score = sil[, 3]
    sil_df = data.frame(cluster=sil[, 1],
                        neighbor=sil[, 2],
                        sil_width=sil[, 3])
    sil_scores_list[[i]] = sil_df
    names(sil_scores_list)[[i]] = paste("res",
                                        as.character(clustering_res[i]),
                                        sep = "_")
  }
  
  # Put sil scores for all resolutions into a single df
  resolutions_sil_df = data.table::rbindlist(sil_scores_list, 
                                             idcol = TRUE)
  return(resolutions_sil_df)
}


#' Internal helper function to get the max sil score from calc_sil_scores()
#' 
#' @param sil_scores_df A data.frame output by calc_sil_scores()
#' @return A numeric indicating the maximum average sil sclore value across tested resolutions
max_avg_sil_score = function(sil_scores_df) {
  avg_sil_scores = sil_scores_df %>%
    group_by(.id) %>%
    summarize(mean_sil = mean(sil_width))
  max_sil_res = as.character(avg_sil_scores[which.max(avg_sil_scores$mean_sil), 1])
  max_sil_res = as.numeric(str_split(max_sil_res, "_")[[1]][2])
  return(max_sil_res)
}


#' Correlation of expression for a target gene across the entire transcriptome. 
#'
#' This function will calculate the correlation between a target gene and other genes in a specific group of cells. 
#' 
#' @param cell_types A character vector indicating cell names of interest.
#' @param metadata_df A data.frame with the metadata for the counts matrix.
#' @param exp_matrix A genes x cells counts matrix.
#' @param target_gene A string indicating the target gene to correlate with the rest of the transcriptome. 
#' @param cor_method A string indicating the correlation method to use, either spearman or pearson
#' 
#' @return A data.frame with the correlation of the target gene with all other genes in the count matrix. 
#' @export
#' 
calc_gene_cors = function(cell_types, metadata_df, 
                          exp_matrix, target_gene, 
                          cor_method) {
  # Create a vector with the names of the cells of interest
  cell_types = cell_types
  metadata = metadata_df
  cells_vector = metadata %>%
    dplyr::filter(prelim_annotations %in% cell_types) %>%
    rownames()
  
  # Process matrix to run correlations
  gene_expression = exp_matrix[, cells_vector]
  reshaped_matrix = t(as.matrix(gene_expression))
  target_gene_vector = as.numeric(reshaped_matrix[, target_gene])
  cors_list = apply(reshaped_matrix, 2, 
                    function(x){cor.test(target_gene_vector, x,
                                         method = cor_method)})
  cors_list_tidy = lapply(cors_list, 
                          broom::tidy)
  cors_df =  data.table::rbindlist(cors_list_tidy)
  gene_names = names(cors_list)
  cors_df_genes = cors_df %>%
    mutate(gene=gene_names,
           FDR = p.adjust(p.value, "BH")) %>%
    arrange(desc(estimate)) %>%
    drop_na() %>%
    mutate(gene_index = seq_along(gene)) %>%
    filter(!gene==target_gene) %>%
    dplyr::rename(cor_estimate = estimate)
  return(cors_df_genes)
}










