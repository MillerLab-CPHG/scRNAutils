# library(Seurat)
# library(tidyverse)
# library(biomaRt)
# library(RColorBrewer)
# library(GOSemSim)
# library(ggsci)



#' Human and mouse homologs
#' 
#' This is a handy function to gene human homologs for mice genes.
#' 
#' @param gene_vec A vector of mice genes
#' @param mouse_biomart A mouse biomart db from which gene attributes will be extracted
#' @return A List where the first item is a data.frame object with 
#' mice gene names, corresponding human homologs, chromosome and type of orthology.
#' Second obj in the list is a vector with the names of genes that have no human homolog and
#' thus could not be converted.
#' @export
convert_mouse_to_human = function(gene_vec, 
                                  mouse_biomart) {
  
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


#' This function will convert gene symbols to Ensembl IDs
#' 
#' @param gene_symbols A character vector with gene symbols to convert to Ensembl IDs. 
#' @param human_biomart A human biomart object.
#' @return A data.frame with gene symbols and their respective Ensembl IDs.
#' 
#' @export
#' 
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


#' This function will process a mice cell type markers df to to keep only genes with one-to-one human orthology and 
#' will add the corresponding human homolog to each row. This will add human homologs in order of descending log2FC. 
#' 
#' @param cell_type_df A data.frame with a list of markers for a cell type. The column name for the gene symbols has to be "Gene_symbol".
#' @param human_homologs_df A data.frame of mouse-human homologs output by covert_mouse_to_human() 
#' @return A data.frame with one-to-one orthologs and the corresponding human homologs.
add_human_homologs_to_mice_df = function(cell_type_df, human_homologs_df) {
  cell_type_df_filtered = cell_type_df[cell_type_df$Gene_symbol %in% human_homologs_df$external_gene_name,]
  idx = match(cell_type_df_filtered$Gene_symbol, human_homologs_df$external_gene_name)
  cell_type_df_filtered$human_homolog = human_homologs_df$hsapiens_homolog_associated_gene_name[idx]
  return(cell_type_df_filtered)
}


#' This is a function to plot gene set enrichment analyses from gProfiler2
#' gprofiler_out: A list containing gProfiler2 output grom gost function
#' @param min_n_genes A numeric indicating Min number of genes in gene set
#' @param max_n_genes A numeric indicating max number of genes in gene set
#' @param pathway_db A character indicating the name of the pathway db to be tested for enrichment
#' @param n_hits A numeric indicating how many terms should be plotted. Terms are ordered in descending FDR
#' @return A ggplot object
#' 
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


#' This function will take a df output with GO terms by gProfiler enrichment analysis
#' and remove GO redundant terms to get a clearer picture of relevant biological processes.
#' @param gprofiler_res_df A data.frame output by gProfiler gost() function. This should be the second item in the output list. 
#' @param cutoff A numeric indicating the cutoff for semantic similarity
#' @param select_fun A function to select one term to keep from a subset of redundant terms based on p values. 
#' @param measure A character indicating the semantic simlarity metric to be used. 
#' @param ontology A character indicating which GO ontology to use (e.g., "BP")
#' @param semData The output of GoSemSim goData function that will allow comparison of semantic similarity between 2 terms. 
#' @return A vector with redundant GO terms meeting the similarity cutoff that should be removed. 
#' 
#' @export 
#' 
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


















