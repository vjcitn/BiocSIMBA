
#' assemble some constants for a run
#' @export
simba_config = function(min_n_genes = NULL, norm_method = 'lib_size', n_top_genes = NULL,
   disc_n_bins = 5L,
   gen_graph_copy = FALSE, gen_graph_hvg = FALSE, gen_graph_dirname = "graph0",
   train_auto_wd = TRUE, train_save_wd = TRUE, train_output = "model",
   train_wd_interval = 10L, train_workers = 12L, train_wd = .015521) {
 list(min_n_genes = min_n_genes, norm_method = norm_method, n_top_genes = n_top_genes, disc_n_bins = disc_n_bins,
   gen_graph_copy = gen_graph_copy, gen_graph_hvg = gen_graph_hvg, gen_graph_dirname = gen_graph_dirname,
   train_auto_wd = train_auto_wd, train_save_wd = train_save_wd, train_output = train_output,
   train_wd_interval = train_wd_interval, train_workers = train_workers, train_wd = train_wd)
}
