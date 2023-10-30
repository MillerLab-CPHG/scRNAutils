



#' Named list of seurat gene expression plots
#' 
#' A function that takes a seurat object as input and 
#' produces a named list of gene expression plots. 
#' 
#' @param seuratObj A seurat object with normalized expression values.
#' @param targetAssay A character indicating the assay to be used for making gene expression plots. 
#' Most of the time this argument will be either RNA or SCT. Default is set to SCT. 
#' @param queryGenes A character vector with the names of the genes to be included in the list.
#' @param splitByVar A character indicating a metadat variable by which to facet the plot (e.g., sex, disease status).
#' @param ptSize A numeric indicating the size of the points for each plot. Default value is 0.1 but can be changed. 
#' 
#' @return A List of named ggplot objects, one for each of the genes included in the character vector.
#' @export
#' 
plotGeneUMAP = function(
  seuratObj,
  queryGenes,
  targetAssay = "SCT", 
  splitByVar = NULL,
  ptSize = 0.1
  ){
  
  # Check original default assay
  originalDefAssay = Seurat::DefaultAssay(seuratObj)
  
  # Set seurat assay to whatever user defines as target assay.
  # If target different than original default assay, do a temporary change. 
  if (originalDefAssay == "RNA" & targetAssay == "SCT") {
    DefaultAssay(seuratObj) = "SCT"
  } else if (originalDefAssay == "SCT" & targetAssay == "RNA") {
    DefaultAssay(seuratObj) = "RNA"
  } else {
    DefaultAssay(seuratObj) = originalDefAssay
  }
  
  plot_list = list()
  if (!is.null(splitByVar)) { 
    gene_plots = lapply(
      queryGenes,
      function(x){FeaturePlot(seuratObj, 
                              features = x, 
                              pt.size = ptSize,
                              split.by = splitByVar,
                              raster = FALSE, 
                              order = TRUE) & custom_theme() & miller_continous_scale()})
  } else {
    gene_plots = lapply(
      queryGenes,
      function(x){FeaturePlot(seuratObj, 
                              features = x, 
                              pt.size = ptSize,
                              raster=FALSE, 
                              order=TRUE) & custom_theme() & miller_continous_scale()})
  }
  
  # We need to return seurat object to the original default assay, 
  # if the target is different of course. 
  if (originalDefAssay != targetAssay) { 
    DefaultAssay(seuratObj) = originalDefAssay
  }
  names(gene_plots) = queryGenes
  return(gene_plots)
}

#' Make a function to visualize the correlation between two genes. 
#' 
#' This function will take a seurat object as input and produce a dot plot 
#' correlating 2 genes in a specific group of cell types.
#' 
#' @param seurat_obj A seurat object with normalized expression values
#' @param target_coldata A character indicating the name of the metadata column with 
#' the cell types of interest.
#' @param target_anno A character vector with the names of the cell types in which we wish 
#' to correlate the genes. 
#' @param assay A character indicating the assay from which to extract normalized counts. One of "SCT" or "RNA".  
#' @param show_cor A boolean indicating whether we wish to calculate and show the correlation coefficient.
#' @param cor_method A character indicating the correlation method to use. One of "pearson" or "spearman". 
#' @return A ggplot object correlating expression of two genes in specific cell populations.  
#'
#' @export 
#' 
plot2GeneCor = function(seurat_obj, 
                        target_coldata,
                        target_anno, target_genes, 
                        show_cor=FALSE, cor_method="pearson",
                        assay="SCT") {
  
  # Tidy df before plotting
  df = scRNAutils::tidy_seurat_counts(seurat_obj=seurat_obj,
                                      target_coldata=target_coldata,
                                      target_anno=target_anno,
                                      assay=assay)
  p = df %>% 
    ggplot(aes_string(x=target_genes[1], 
                      y=target_genes[2])) + 
    geom_point() + 
    custom_theme()  
  
  # Show correlation coefficient if the user wants to do so
  if(show_cor) { 
    # Fit linear model 
    p = p + 
      geom_smooth(method="lm") + 
      ggpubr::stat_cor(method=cor_method)
  }
  
  return(p)
}

#' Make a function to plot correlations
#' 
#' This function will take the output of calc_gene_cors() and produce a
#' nice way to visualize the correlations. 
#' 
#' @param cors_df A data.frame that is output by calc_gene_cors()
#' @param target_gene A character string indicating the target gene
#' @param target_cells A character indicating the cells where the cors were calculated 
#' @param genes_to_label A character vector with the names of the genes that should be highlighted in the plot.
#' @param point_size A numeric indicating the size of the points in the plot  
#' @return A ggplot object
#' 
#' @export
#' 
plotGene2TranscriptomeCor = function(
  cors_df, 
  target_gene, 
  target_cells, 
  genes_to_label, 
  point_size=0.9
  ){
  
  ggplot(cors_df, aes(x=gene_index, y=cor_estimate, 
                      color=cor_estimate, 
                      label=ifelse(gene %in% genes_to_label, gene, ""))) + 
    geom_point(size=point_size) +
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
    scale_x_continuous(limits=c(0, 20000)) +
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(face="italic", size=20))
  #legend.position = "bottom")
}


#' Model gene expression across pseudotime. 
#' 
#' This function will plot normalized expression values across pseudotime
#' 
#' @param cds_ordered: A monocle3 cds object with cells ordered across pseudotime.
#' @param cell_annotations: A vector with cell types of interest for plotting gene expression.
#' @param genes A Vector with the genes to be plotted across pseudotime.
#' @param facet_wrap_plot A boolean indicating whether all genes shoud be shown in one or multiple plots.
#' @param pseudotime_boundary A numeric indicating the upper boundary of pseudotime to be plotted. 
#'  
#' @return A ggplot object with one or multiple panels showing gene expression across pseudotime values.
#' @export
#' 
plot_expression_on_pseudotime = function(
  cds_ordered, 
  cell_annotations, 
  genes,
  facet_wrap_plot = FALSE,
  pseudotime_boundary = 30
  ){ 
  
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
    geom_smooth(method = "lm", 
                formula = y ~ splines::ns(x, 3), 
                color="black") 
  p = p + 
    xlim(0, pseudotime_boundary) + 
    scale_y_log10() + 
    ylab("SCTransform norm expression") + 
    custom_theme()
  
  # If we want to plot genes on different panels
  if (facet_wrap_plot) { 
    p = p + 
      facet_wrap(~Feature_id, scales = "free_y") + 
      theme(legend.position = "none",
            strip.background = element_rect(fill="white"))
    return(p)
  }
  return(p)
}

#' Internal helper function to plot scRNA library complexity
#'
#' this helper function will correlate number of UMIs to
#' the number of genes detected for each cell, and color them
#' by the percent of reads aligned to the mitochondrial genome
#'
#' @param seuratObj A Seurat object with QC metrics stored in the metadata 
#' @param df A dataframe with metadata variables nCount_RNA, nFeature_RNA
#' and percent.mt variables.
#' @param cor A boolean indicating whether we want to apply a lm to 
#' UMIs and genes. 
#' @return A ggplot object showing the correlation between 
#' number of UMIs and detected genes. 
#'
#' @export
#'
plotRNAcomplexity = function(
  seuratObj,
  df,
  cor = FALSE
  ){
  df = seuratObj@meta.data
  p = ggplot(df, 
             aes_string(x = "nCount_RNA", y = "nFeature_RNA", 
                 color = "percent.mt"))
  
  p = p +
    geom_point(size = 0.4) +
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("nUMIs") +
    ylab("nGenes") + 
    ggtitle("Sample complexity") + 
    custom_theme() +
    miller_continous_scale()
  
  if (cor) { 
    p = p + 
      geom_smooth(method = "lm") + 
      stat_cor()
  }
  p 
}

#' Helper function to plot features
#' 
#' This function is equivalent to Seurat vlnPlot
#' but can be customized in any way we want to show
#' important QC metrics such as number of UMIs, detected
#' genes or percentage of reads mapped to mitochondrial 
#' or ribsomal genes.
#' 
#' @param df A dataframe containing information about QC metrics.
#' This df should have a "libraryID" column.
#' @param feature A character vector indicating the metadata variable 
#' to plot (i.e. nUMIs)
#' 
#' @return A ggplot object showing the distribution of QC metrics
#' across each cell in the library. 
#' 
#' @export
#' 
plotRNAFeature = function(seuratObj, feature) { 
  df = seuratObj@meta.data
  ggplot(df, aes_string(x="sample", y=feature)) + 
    geom_violin(width = 0.3, fill = "azure3") + 
    #geom_boxplot(width=0.2, outlier.shape = NA) + 
    geom_jitter(width = 0.04, alpha=0.2, size=0.4) + 
    custom_theme() +   
    theme(aspect.ratio = 1.3)
  
}

#' Plot features list
#' 
#' Wrapper of `plotRNAFeatures` to show several features
#' at the same time
#' 
#' @param seuratObj A Seurat object with metadata relevant to
#' QC metrics.
#' @param features A character vector with several QC metrics
#' to plot. 
#' 
#' @xport A list of ggplot objects with any specified metrics. 
#' 
plotRNAFeatureList = function(seuratObj, features) { 
  features = features
  plotList = lapply(features, plotRNAFeature, 
                    seuratObj = seuratObj)
  metricsPlot = ggarrange(plotlist = plotList)
  return(metricsPlot)
  }

#' Custom theme to make ggplots prettier
#' 
#' I'm not a big fan of the standard way Seurat ouputs UMAP plots
#' so this function will provide a proper aspect ratio and also
#' decently sized x and y axis label. 
#' 
#' @return A ggplot theme 
#' @export
custom_theme =  function() {
  return(theme_bw() + 
           theme(
             aspect.ratio = 1,
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             axis.text = element_text(size=12),
             axis.title = element_text(size=12),
             title = element_text(size=12),
             legend.text = element_text(size=12)
           )
  )
}


#' A nicer continuous scale than the default ggplot one. 
#'
#'@param style A character string indicating whether we need to plot points or bars. 
#'
#'@return A ggplot continuous scale
#'@export
#'
miller_continous_scale = function(style="points") {
  
  color_palette = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF",
                    "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#053061")
  
  if (style=="bars") {
    cont_scale = scale_fill_gradientn(colours = rev(color_palette))
  } else {
    cont_scale = scale_colour_gradientn(colours = rev(color_palette))
  }
  return(cont_scale)
}


#' A nicer discrete scale than the ggplot default one. 
#'
#' @param style A character string indicatin whether to plot points or bars.
#' @param option A numeric indicating whether to plot the first or second option.
#' 
#' @return A ggplot discrete scale
#' 
#' @export  
#'
miller_discrete_scale = function(style="points", option=1) { 
  npg = pal_npg("nrc")(10)
  nature_scale = c(npg, "darkgoldenrod1")
  nature_scale2 = c("#6A3D9A", nature_scale[-1])
  
  if (style=="bars") {
    if (option==1) {
      discr_scale = scale_fill_manual(values = nature_scale)
    } else if (option==2) {
      discr_scale = scale_fill_manual(values = nature_scale2)
    }
  } else {
    if (option==1) {
      discr_scale = scale_colour_manual(values = nature_scale)
    } else if (option==2) {
      discr_scale = scale_colour_manual(values = nature_scale2)
    }
  }

  return(discr_scale)
}



