#simba_config = function(min_n_genes = NULL, norm_method = 'lib_size', n_top_genes = NULL,
#   disc_n_bins = 5L,
#   gen_graph_copy = FALSE, gen_graph_hvg = FALSE, gen_graph_dirname = "graph0",
#   train_auto_wd = TRUE, train_save_wd = TRUE, train_output = "model",
#   train_wd_interval = 10L, train_workers = 12L, train_wd = .015521) {


#' train Simba PBG embedding for an SCE (or h5ad)
#' @examples
#' h5adpath = "/home/vincent/tenx3k.h5ad"
#' train_sce(h5adpath)
#' @export
train_sce = function(h5adpath, simconf = simba_config()) {
  proc = basilisk::basiliskStart(bsklenv, testload="simba") # avoid package-specific import
  basilisk::basiliskRun(proc, function(h5ad, simba_config) {
     sref = reticulate::import("simba")
     adata_CG = sref$read_h5ad(h5adpath)
     conf = simconf
     if (!is.null(conf$min_n_genes)) sref$pp$filter_cells_rna(adata_CG,
               as.integer(conf$min_n_genes))
     sref$pp$normalize(adata_CG, method = conf$norm_method)
     sref$pp$log_transform(adata_CG)
     if (!is.null(conf$n_top_genes)) sref$pp$select_variable_genes(adata_CG,
               as.integer(conf$n_top_genes))
     sref$tl$discretize(adata_CG, n_bins = as.integer(conf$disc_n_bins))
    # generate graph
     sref$tl$gen_graph(list_CG=list(adata_CG), copy=conf$gen_graph_copy,
	use_highly_variable = conf$gen_graph_hvg, dirname = conf$gen_graph_dirname) 
     sref
     }, h5adpath, simconf)
}

