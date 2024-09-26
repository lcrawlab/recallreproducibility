#' @title Given a seurat object, computes standard clustering metrics for clustering assignments vs true group labels.
#'
#' @description Given a seurat object, computes standard clustering metrics for clustering assignments vs true group labels.
#'
#' @param seurat_obj The Seurat object containing clustering results.
#' @param truth_col The name of the metadata column containing the true groups.
#' @param cluster_label_name The name of the metadata column containing the clusters.
#' @param replicate The replicate number (for storing in the output dataframe).
#' @param num_cells The number of cells (for storing in the output dataframe).
#' @param num_groups The number of groups (for storing in the output dataframe).
#'
#' @returns A dataframe with the metrics.
#'
#' @name get_clustering_metrics
#' @export
get_clustering_metrics <- function(seurat_obj,
                                   truth_col,
                                   cluster_label_name,
                                   replicate,
                                   num_cells,
                                   num_groups) {

  ari <- pdfCluster::adj.rand.index(seurat_obj@meta.data[[truth_col]], 
                                    seurat_obj@meta.data[[cluster_label_name]])

  v_measure <- clevr::v_measure(seurat_obj@meta.data[[truth_col]], 
                                        seurat_obj@meta.data[[cluster_label_name]])


  homogeneity <- clevr::homogeneity(seurat_obj@meta.data[[truth_col]], 
                                            seurat_obj@meta.data[[cluster_label_name]])
  
  completeness <- clevr::completeness(seurat_obj@meta.data[[truth_col]], 
                                              seurat_obj@meta.data[[cluster_label_name]])
  
  fowlkes_mallows <- clevr::fowlkes_mallows(seurat_obj@meta.data[[truth_col]], 
                                                    seurat_obj@meta.data[[cluster_label_name]])
  jaccard <- clusteval::jaccard(seurat_obj@meta.data[[truth_col]], 
                                        seurat_obj@meta.data[[cluster_label_name]])
  
  num_clusters <-  length(levels(as.factor(seurat_obj@meta.data[[cluster_label_name]])))
  method <- cluster_label_name

  
  metrics_df_row <- data.frame(method,
                               replicate,
                               num_cells,
                               num_groups,
                               num_clusters,
                               ari,
                               v_measure,
                               fowlkes_mallows,
                               homogeneity,
                               completeness,
                               jaccard)
  return(metrics_df_row)
}
