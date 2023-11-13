#simba_config = function(min_n_genes = NULL, norm_method = 'lib_size', n_top_genes = NULL,
#   disc_n_bins = 5L,
#   gen_graph_copy = FALSE, gen_graph_hvg = FALSE, gen_graph_dirname = "graph0",
#   train_auto_wd = TRUE, train_save_wd = TRUE, train_output = "model",
#   train_wd_interval = 10L, train_workers = 12L, train_wd = .015521) {


#' produce graph underlying Simba PBG embedding for an sc-RNA-seq SCE (or h5ad)
#' @param h5adpath character(1) path to h5ad
#' @param simconf output of `simba_config`
#' @return an instance of S3 class `simba_scrna`,
#' a list with elements gdf (a data.frame of source, relation, destination, only
#' produced if gen_graph_copy is TRUE in simba_config), work_contents (full pathnames
#' of pbg outputs), workdir (a copy of the simconf$workdir),
#' C_emb (a SingleCellExperiment instance with cell embedding),
#' G_emb (a SingleCellExperiment instance with gene embedding),
#' pp_CG (a SingleCellExperiment instance with preprocessed representation
#' of the input data).
#' @note simba module is imported with convert=FALSE to help zellkonverter.  In
#' the example gout$gdf is a python reference that could go stale.
#' @examples
#' h5adpath = "/home/vincent/tenx3k.h5ad"
#' gout = gen_graph_sce(h5adpath, simconf = simba_config(gen_graph_copy=TRUE)) # generally want copy to be FALSE
#' gout$gdf
#' dir(gout$work_contents, full.names=TRUE, recursive=TRUE)
#' gout
#' @export
gen_graph_sce = function (h5adpath, simconf = simba_config(gen_graph_copy = FALSE)) 
{
    proc = basilisk::basiliskStart(bsklenv, testload = "simba")
    basilisk::basiliskRun(proc, function(h5ad, simba_config) {
        conf = simconf
        sref = reticulate::import("simba", convert=FALSE)
        sref$settings$set_workdir(conf$workdir)
        adata_CG = sref$read_h5ad(h5adpath)
        if (!is.null(conf$min_n_genes)) 
            sref$pp$filter_cells_rna(adata_CG, as.integer(conf$min_n_genes))
        sref$pp$normalize(adata_CG, method = conf$norm_method)
        sref$pp$log_transform(adata_CG)
        if (!is.null(conf$n_top_genes)) 
            sref$pp$select_variable_genes(adata_CG, as.integer(conf$n_top_genes))
        sref$tl$discretize(adata_CG, n_bins = as.integer(conf$disc_n_bins))
        gdf = sref$tl$gen_graph(list_CG = list(adata_CG), copy = conf$gen_graph_copy, 
            use_highly_variable = conf$gen_graph_hvg, dirname = conf$gen_graph_dirname)
        sref$tl$pbg_train(auto_wd = TRUE, save_wd = TRUE, output = "model")
        ndict = sref$read_embedding()
        C_emb = zellkonverter::AnnData2SCE(ndict["C"], hdf5_backed=FALSE)
        G_emb = zellkonverter::AnnData2SCE(ndict["G"], hdf5_backed=FALSE)
        origAD = zellkonverter::AnnData2SCE(adata_CG)
        work_contents = dir(conf$workdir, full.names = TRUE)
        ans = list(gdf = gdf, work_contents = work_contents, workdir = conf$workdir, 
            C_emb = C_emb, G_emb = G_emb, ppCG=origAD)
        class(ans) = c("simba_scrna", "list")
        ans
    }, h5adpath, simconf)
}

#' print method
#' @param x instance of simba_scrna
#' @param \dots not used
#' @export
print.simba_scrna = function(x, ...) {
  dd = dim(x$ppCG)
  cat(sprintf("simba_scrna instance based on input with %d genes and %d cells.\n",
        dd[1], dd[2]))
  cat(sprintf("  The workdir used was %s.\n", x$workdir))
}


gen_graph_sce_bad = function(h5adpath, simconf = simba_config(gen_graph_copy=FALSE)) {
  proc = basilisk::basiliskStart(bsklenv, testload="simba") # avoid package-specific import
  basilisk::basiliskRun(proc, function(h5ad, simba_config) {
     conf = simconf
     sref = reticulate::import("simba")
     sref$settings$set_workdir(conf$workdir)
     adata_CG = sref$read_h5ad(h5adpath)
     if (!is.null(conf$min_n_genes)) sref$pp$filter_cells_rna(adata_CG,
               as.integer(conf$min_n_genes))
     sref$pp$normalize(adata_CG, method = conf$norm_method)
     sref$pp$log_transform(adata_CG)
     if (!is.null(conf$n_top_genes)) sref$pp$select_variable_genes(adata_CG,
               as.integer(conf$n_top_genes))
     sref$tl$discretize(adata_CG, n_bins = as.integer(conf$disc_n_bins))
    # generate graph
     gdf = sref$tl$gen_graph(list_CG=list(adata_CG), copy=conf$gen_graph_copy,
	use_highly_variable = conf$gen_graph_hvg, dirname = conf$gen_graph_dirname) 
     sref$tl$pbg_train(auto_wd = TRUE, save_wd = TRUE, output = 'model' ) # produces files under 'model'
     dict = sref$read_embedding()
     C_emb = zellkonverter::AnnData2SCE(dict['C'])
     G_emb = zellkonverter::AnnData2SCE(dict['G'])
     work_contents = dir(conf$workdir, full.names=TRUE)
     list(gdf = gdf, work_contents=work_contents, workdir=conf$workdir, C_emb=C_emb, G_emb=G_emb)
     }, h5adpath, simconf)
}

