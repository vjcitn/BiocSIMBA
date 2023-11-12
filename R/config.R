
#' assemble some constants for a run
#' @export
simba_config = function(workdir=paste0(tempfile(), "/simba_work"), min_n_genes = NULL, norm_method = 'lib_size', n_top_genes = NULL,
   disc_n_bins = 5L,
   gen_graph_copy = FALSE, gen_graph_hvg = FALSE, gen_graph_dirname = "graph0") {
 list(workdir=workdir, min_n_genes = min_n_genes, norm_method = norm_method, n_top_genes = n_top_genes, disc_n_bins = disc_n_bins,
   gen_graph_copy = gen_graph_copy, gen_graph_hvg = gen_graph_hvg, gen_graph_dirname = gen_graph_dirname)
}



