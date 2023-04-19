#library(Seurat)
#library(tidyverse)
#library(biomaRt)
#library(RColorBrewer)
#library(GOSemSim)
#library(ggsci)

# This script will contain several functions to stremaline scRNA datasets processing
# including removal of doublets/ambient mRNA, SCT normalization and UMAP embeddings. 

# Genes for regressing out cell cycle variance
cc.genes.updated.2019
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

###########################
# Continuous color scales #
###########################

# Options for color scales: spectral, viridis(magma, inferno)
# A new scale for gene expression plots
new_scale = scale_colour_gradientn(colours = rev(brewer.pal(n = 11, 
                                                            name = "Spectral")))
scale_inferno = scale_color_viridis_c(option="inferno")

new_scale_bars = scale_fill_gradientn(colours = rev(brewer.pal(n=11,
                                                          name = "Spectral")))

# Define continous scale for expression UMAPs
spectral_scale = brewer.pal(11, "Spectral")
#spectral_scale = c(spectral_scale, "black") #053061
spectral_scale2 = c(spectral_scale[1:10], "#053061")
new_scale3 = scale_colour_gradientn(colours =  rev(spectral_scale2))
new_scale3_bars = scale_fill_gradientn(colours =  rev(spectral_scale2))

# Color scale for heatmaps
palette_length = 100
blue_red_palette = colorRampPalette(c("#053061", "white", "darkred"))(palette_length)
blue_red_scale = scale_colour_gradientn(colours = blue_red_palette)

#########################
# Discrete color scales #
#########################

# Add other scales to make plots look prettier
# Import color scales
npg = pal_npg("nrc")(10)
aaas = pal_aaas()(10)
nature_scale = c(npg, "darkgoldenrod1")
npg_scale_bars = scale_fill_manual(values=nature_scale)
npg_scale = scale_colour_manual(values = nature_scale)
aaas_scale = c(aaas, "darkgoldenrod1")

# New npg scale
nature_scale2 = c("#6A3D9A", nature_scale[-1])
#nature_scale2[10] = "#E7298A"
npg_scale2 = scale_colour_manual(values = nature_scale2)
npg_scale2_bars = scale_fill_manual(values = nature_scale2)

# Add color scale for level2 annotations
paired = brewer.pal(n=12, name="Paired")
paired = paired[-6]
set2 = brewer.pal(n=8, name="Set2")
dark2 = brewer.pal(n=8, name="Dark2")
large_discrete_scale = c(nature_scale2, set2[1:6], dark2[1:7], paired[1:8])
large_discrete_scale = rev(large_discrete_scale)
large_discrete_scale[22] = "#E7298A"
large_discrete_scale[12] = "#980043"
large_discrete_scale[23] = "#AE017E"

level2_annotations_scale = scale_colour_manual(values=large_discrete_scale)
level2_annotations_scale_bars = scale_fill_manual(values=large_discrete_scale)

# Interpolate nature scale
new_length = 32
nature_extended = colorRampPalette(nature_scale2)(new_length)
nature_ext_scale = scale_colour_manual(values = rev(nature_extended))


#' Internal helper function to automate doublet removal workflow and update metadata
#' @param  Arg seurat_obj: pre-processed seurat object with cell clusters annotations
#' @param nrep: numeric; number of times to run the scDblFinder function
#' @param multisample_dataset: Logical indicating whether the processed dataset has multiple samples
#' @param sample_ids if multisample_dataset = TRUE, A character that refers to the metadata column
#' that contains sample ids
#' @return A vector of consensus doublet cell barcodes
scDblFinder_clusters = function(seurat_obj, nrep, 
                                multisample_dataset=FALSE, sample_ids=NULL) {
  
  sample_sce = as.SingleCellExperiment(seurat_obj)
  
  # If we have more than one sample in a dataset, use the samples parameter
  if (multisample_dataset) {
    sample_sce_dbl = replicate(nrep, scDblFinder(sample_sce,
                                                 samples = sample_ids,
                                                 clusters = "seurat_clusters"))
  } else {
    sample_sce_dbl = replicate(nrep, scDblFinder(sample_sce, 
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


#' Internal helper function to standardize removal of ambient mRNA contamination
#' @param seurat_obj: A seurat obj that has been filtered for doublets
#' @return A seurat object with a decontaminated raw counts matrix

decontX_remove = function(seurat_obj) {
  decontX_results = celda::decontX(seurat_obj@assays$RNA@counts)
  decontaminated_matrix = decontX_results$decontXcounts
  seurat_obj@assays$RNA@counts = decontaminated_matrix
  return(seurat_obj)
}


#' Internal helper function to standardize filtering and SCT normalization. We might need 
#' to leave the FindClusters() function separate since that really depends on the
#' number of cells on each dataset. 
#' @param seurat_obj An input seurat object.
#' @param seurat_filter A boolean indicating whether do we wish to filter. 
#' Set to FALSE in pre-processing round for doublet removal
#' @param sample A character indicating name of the sample.
#' @param study A character indicating name of the study. 
#' @return A seurat object that has been SCT normalized, nearest neighbors graph and UMAP embeddings
Seurat_SCT_process = function(seurat_obj, seurat_filter=FALSE, 
                              sample_id, study_name, artery, disease_status){
  
  # Define "stim" group to separate datasets during visualization after integration. 
  # It might be helpful to plot by sample ID, study. Vascular bed and sample disease status will be added once we create the main reference since it's a bit trickier to set
  seurat_obj$sample = sample_id
  seurat_obj$study = study_name
  
  
  # Check for mt, hb percentage and other quality metrics 
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Define hemoglobin genes
  hb_index = grep(rownames(seurat_obj), pattern = "^HB[AB]")
  hb_genes = rownames(seurat_obj)[hb_index]
  seurat_obj[["percent.hb"]] = PercentageFeatureSet(seurat_obj, features = hb_genes)
  
  # Filter low quality cells (start with parameters used in the paper)
  # NOTE: Skip this step during the first processing round to be used for doublet detection
  if (seurat_filter) {
    seurat_obj =  subset(seurat_obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & nCount_RNA >= 200 & nCount_RNA <= 20000 & percent.mt <= 10 & percent.hb <= 5)
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



#' Helper function to calculate clusters silhouette scores from a seurat object
#' @param seurat_obj A seurat obj with reduced dims embedding like pca. Default now is set to PCA. 
#' @param clustering_res A numeric vector with clustering resolutions to test
#' @return A data.frame with silhouette scores for a the provided clustering resolutions
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
  dist_matrix = dist(x = Seurat::Embeddings(object = seurat_obj[["pca"]])[, 1:30])
  
  # Cluster data with each of the resolutions provided
  for (i in seq_along(clustering_res)) { 
    seurat_obj = FindClusters(seurat_obj, resolution = clustering_res[i])
    
  clusters = seurat_obj$seurat_clusters
  
  message("Calculating silhouette scores for defined resolutions...")
  sil = silhouette(x = as.numeric(x = as.factor(x = clusters)), 
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



# This is a handy function to gene human homologs for mice genes
# Arg gene_vec: A vector of mice genes
# Arg mouse_biomart: mouse biomart db from which gene attrbutes will be extracted
# Value A list where the first item is a data.frame object with 
# mice gene names, corresponding human homologs, chromosome and type of orthology.
# Second obj in the list is a vector with the names of genes that have no human homolog and
# thus could not be converted.
convert_mouse_to_human = function(gene_vec, mouse_biomart) {
  
  # Get attributes of interest from mouse biomart object
  mouse_attributes = c("external_gene_name",
                       "hsapiens_homolog_associated_gene_name", 
                       "hsapiens_homolog_chromosome",
                       "hsapiens_homolog_orthology_type",
                       "hsapiens_homolog_perc_id")
  
  # Get human mouse orthologs
  # This is getting all of the orthologs. First filter for genes that have mouse
  # orthologs
  mouse_human_orthologs = getBM(attributes = mouse_attributes,
                                filters = "with_hsapiens_homolog",
                                values = TRUE,
                                mart = mouse_biomart,
                                uniqueRows = TRUE)
  
  # Filter for genes that have one-to-one human-mouse orthologs
  filtered_orth = mouse_human_orthologs %>%
    filter(hsapiens_homolog_orthology_type == "ortholog_one2one")
  
  # Select user genes that are one-to-one orthologs
  ol_idx = which(gene_vec %in% filtered_orth$external_gene_name)
  user_one_to_one = gene_vec[ol_idx]
  
  # Would be interesting know how many couldn't get converted
  unconverted_genes = setdiff(gene_vec, user_one_to_one)
  
  # Filter one-to-one ortholog list
  filtered_orth = filtered_orth %>%
    filter(external_gene_name %in% user_one_to_one)
  
  return(list(filtered_orth, unconverted_genes))
}



# This function will conver gene symbols to Ensembl IDs
# Args gene_symbols: A vector with gene symbols to convert to Ensembl IDs 
# Args human_biomart: A human biomart object
# Value: A df with gene symbols and their respective Ensembl IDs

get_human_ensembl_ids = function(gene_symbols, human_biomart) {
  
  human_attributes = c("external_gene_name",
                       "ensembl_gene_id")
  
  ensembl_ids = getBM(attributes = human_attributes,
                      values = TRUE,
                      mart = human_biomart,
                      uniqueRows = TRUE)
  
  target_ens_ids = ensembl_ids[ensembl_ids$external_gene_name %in% gene_symbols, ]
  
  # There are several Ensembl IDs for some genes so we'll remove duplicates 
  target_ens_ids_no_dups = target_ens_ids[!duplicated(target_ens_ids$external_gene_name), ]
  
  return(target_ens_ids_no_dups)
  
}







# This function will process a mice cell type markers df to to keep only genes with one-to-one human orthology and 
# will add the corresponding human homolog to each row. This will add human homologs in order of descending log2FC. 
# Args cell_type_df: a data.frame with a list of markers for a cell type. The column name for the gene symbols has to be "Gene_symbol"
# Args human_homologs_df: a data.frame of mouse-human homologs output by covert_mouse_to_human() 

add_human_homologs_to_mice_df = function(cell_type_df, human_homologs_df) {
  cell_type_df_filtered = cell_type_df[cell_type_df$Gene_symbol %in% human_homologs_df$external_gene_name,]
  idx = match(cell_type_df_filtered$Gene_symbol, human_homologs_df$external_gene_name)
  cell_type_df_filtered$human_homolog = human_homologs_df$hsapiens_homolog_associated_gene_name[idx]
  return(cell_type_df_filtered)
}

# This function will calculate the correlation between a target gene and other genes in a specific group of cells. 
# Make a function to calculate correlations between a target gene and genes in a specific group of cells
calc_gene_cors = function(cell_types, metadata_df, exp_matrix, target_gene, cor_method) {
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
                    function(x){cor.test(target_gene_vector, x, method = cor_method)})
  cors_list_tidy = lapply(cors_list, tidy)
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

# Make a function to plot correlations
plot_cors = function(cors_df, target_gene, target_cells, genes_to_label) {
  ggplot(cors_df, aes(x=gene_index, y=cor_estimate, 
                      color=cor_estimate, label=ifelse(gene %in% genes_to_label, gene, ""))) + 
    geom_point(size=0.9) +
    geom_text_repel(min.segment.length = 0, 
                    max.overlaps = Inf, color="black", 
                    force=1, fontface="italic", size=4, 
                    segment.color = "grey27") +
    ggtitle(paste0(target_gene, " Pearson correlations in ", target_cells)) +
    xlab("Genes") + 
    ylab("Cor coeff") + 
    labs(color="r") +
    theme_bw() +
    scale_colour_viridis_c(option="magma") +
    scale_x_continuous(limits=c(0,20000)) +
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(face="italic", size=20))
  #legend.position = "bottom")
}


# This function will plot normalized expression values across pseudotime
# Args cds_ordered: A monocle3 cds object with cells ordered across pseudotime
# Args cell annotations: A vector with cell types of interest for plotting gene expression
# Args genes: A vector with the genes to be plotted across pseudotime
# Args facet_wrap_plot: A boolean indicating whether all genes shoud be shown in one or multiple plots
# Value: A ggplot object with one or multiple panels showing gene expression across pseudotime values. 
plot_expression_on_pseudotime = function(cds_ordered, cell_annotations, genes,
                                         facet_wrap_plot=FALSE,
                                         pseudotime_boundary = 30) { 
  
  # Create subset of the main cds object to contain only cell type annotations
  # and genes for plotting along pseudotime. 
  cds_subset = cds_ordered[rowData(cds_ordered)$gene_short_name %in% genes,
                           colData(cds_ordered)$prelim_annotations %in% cell_annotations]
  
  # Create df of subset for gene expression across pseudotime plotting
  cds_exprs = SingleCellExperiment::counts(cds_subset)
  cds_exprs = Matrix::t(Matrix::t(cds_exprs))
  cds_exprs = reshape2::melt(round(as.matrix(cds_exprs)))
  colnames(cds_exprs) = c("Feature_id", "Cell", "Norm_expression")
  
  # Add pseudotime values for each cell
  pseudotime = pseudotime(cds_subset)
  cds_exprs$pseudotime = pseudotime[match(cds_exprs$Cell, 
                                          names(pseudotime))]
  
  # Add annotations (column name of the annotations is fixed
  # for the smc_rpca_sct_v3 subset)
  annotations = colData(cds_subset)$prelim_annotations
  names(annotations) = rownames(colData(cds_subset))
  cds_exprs$annotations = annotations[match(cds_exprs$Cell, 
                                            names(annotations))]
  
  # Order panels according to order of genes provided
  cds_exprs$Feature_id = factor(cds_exprs$Feature_id,
                                levels = genes)
  
  # Fit a cubic spline to normalized expression values across pseudotime
  p =  ggplot(cds_exprs, aes(pseudotime, Norm_expression, fill=Feature_id)) + 
    geom_smooth(method = "lm", formula = y ~ splines::ns(x, 3), 
                color="black") 
  p = p + 
    xlim(0, pseudotime_boundary) + 
    scale_y_log10() + 
    ylab("SCTransform norm expression") + 
    theme_bw() + 
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    npg_scale_bars
  
  # If we want to plot genes on different panels
  if (facet_wrap_plot) { 
    p = p + 
      facet_wrap(~Feature_id, scales = "free_y") + 
      theme(legend.position = "none",
            strip.background = element_rect(fill="azure3"))
    return(p)
  }
  return(p)
}






##################################################
# gProfiler2 and other plotting helper functions #
##################################################

# This is a function to plot gene set enrichment analyses from gProfiler2
# Arg gprofiler_out: A list containing gProfiler2 output grom gost function
# Arg min_n_genes: A numeric indicating Min number of genes in gene set
# Arg max_n_genes: A numeric indicating max number of genes in gene set
# Arg pathway_db: A character indicating the name of the pathway db to be tested for enrichment
# Arg n_hits: A numeric indicating how many terms should be plotted. Terms are ordered in descending FDR

gprofiler_bar_plot = function(gprofiler_df, min_n_genes, max_n_genes, 
                              pathway_db, n_hits) {
  gprofiler_df %>%
    filter(term_size >= min_n_genes & term_size <= max_n_genes & source == pathway_db) %>%
    head(n=n_hits) %>%
    mutate(log_padj = -log10(p_value)) %>%
    arrange(desc(log_padj)) %>%
    ggplot(aes(x=reorder(term_name, log_padj), y=log_padj)) + 
    geom_col(width = 0.4, fill="#053061") +
    theme_bw() +
    ylab("-log10(FDR)") +
    xlab(paste(pathway_db, "term", sep = " ")) + 
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkred") +
    coord_flip() +
    theme(aspect.ratio = 1,
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          title = element_text(size=14)) + 
    scale_fill_viridis_c()
}


# This function will take a df output with GO terms by gProfiler enrichment analysis
# and remove GO redundant terms to get a clearer picture of relevant biological processes.
# Arg gprofiler_res_df: A data.frame output by gProfiler gost() function. This should be the second item in the output list. 
# Arg cutoff: A numeric indicating the cutoff for semantic similarity
# Arg select_fun: function to select one term to keep from a subset of redundant terms based on p values. 
# Arg measure: A character indicating the semantic simlarity metric to be used. 
# Arg ontology: A character indicating which GO ontology to use (e.g., "BP")
# Arg semData: The output of GoSemSim goData function that will allow comparison of semantic similarity between 2 terms. 
# Value: A vector with redundant GO terms meeting the similarity cutoff that should be removed. 
# Example: 
# hsGO = GOSemSim::godata("org.Hs.eg.db, ont="BP")
# c8_gprofiler_res = c8_gprofiler$result
# c8_go_bp = c8_go_bp[c8_go_bp$source == "BP", ]
# redundant_terms = gprofiler_go_simplify(c8_go_bp, semData=hsGO, ontology= "BP")

gprofiler_go_simplify = function(gprofiler_res_df, cutoff=0.7, select_fun=min, 
                             measure="Rel", ontology, semData) {
  if (missing(semData) || is.null(semData)) {
    if (measure == "Wang") {
      semData = godata(ont = ontology)
    } else {
      stop("godata should be provided for IC-based methods...")
    }
  } else {
    if (ontology != semData@ont) {
      msg = paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
      stop(msg)
    }
  }
  
  sim = mgoSim(gprofiler_res_df$term_id, gprofiler_res_df$term_id,
               semData = semData,
               measure=measure,
               combine=NULL)
  
  ## to satisfy codetools for calling gather
  go1 = go2 = similarity <- NULL
  
  sim.df = as.data.frame(sim)
  sim.df$go1 = row.names(sim.df)
  sim.df = gather(sim.df, go2, similarity, -go1)
  
  sim.df = sim.df[!is.na(sim.df$similarity),]
  
  ## feature 'by' is attached to 'go1'
  sim.df = merge(sim.df, gprofiler_res_df[, c("term_id", "p_value")], by.x="go1", by.y="term_id")
  sim.df$go2 = as.character(sim.df$go2)
  
  ID = gprofiler_res_df$term_id
  
  GO_to_remove = character()
  for (i in seq_along(ID)) {
    ii = which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    ## if length(ii) == 1, then go1 == go2
    if (length(ii) < 2)
      next
    
    sim_subset = sim.df[ii,]
    
    jj = which(sim_subset[, "p_value"] == select_fun(sim_subset[, "p_value"]))
    
    ## sim.df <- sim.df[-ii[-jj]]
    GO_to_remove = c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }
  
  #res[!res$ID %in% GO_to_remove, ]
  return(GO_to_remove)
}


# Write a function that takes a seurat object as input and 
# produces a list of gene expression plots 
# Write a function to generate feature plots and store them within a list
# Args seurat_object A seurat object with normalized expression values
# Args genes_of_interest A character vector with the names of the genes to be included in the list
# Args pt.size A numeric indicating the size of the points for each plot. Default value is 0.1 but can be changed. 
# Value A list of ggplot objects, one for each of the genes included in the character vector
# What if we want to split plots by a variable like disease status
seurat_gene_plot_list = function(seurat_object, assay="SCT", 
                                 genes_of_interest,
                                 split_by_var=NULL,
                                 pt.size=0.1) {
  if (assay=="RNA") {
    DefaultAssay(seurat_object) = "RNA"
  }
  plot_list = list()
  if (!is.null(split_by_var)) { 
    gene_plots = lapply(genes_of_interest,
                        function(x){FeaturePlot(seurat_object, 
                                                features = x, pt.size = pt.size,
                                                split.by = split_by_var,
                                                raster = FALSE, 
                                                order = TRUE)} & custom_theme & new_scale3)
  } else {
    gene_plots = lapply(genes_of_interest,
                      function(x){FeaturePlot(seurat_object, 
                        features = x, pt.size = pt.size,
                        raster = FALSE, 
                        order = TRUE)} & custom_theme & new_scale3)
  }
  names(gene_plots) = genes_of_interest
  return(gene_plots)
}


# This is a custom theme to make ggplots prettier
custom_theme =  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        title = element_text(size=14),
        legend.text = element_text(size=14))











