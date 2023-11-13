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
#' h5adpath = system.file(file.path("h5ad", "tenx3k.h5ad"), package="BiocSIMBA")
#' gout = gen_graph_sce(h5adpath, simconf = simba_config(gen_graph_copy=TRUE)) # generally want copy to be FALSE
#' gout$gdf
#' dir(gout$work_contents, full.names=TRUE, recursive=TRUE)
#' gout
#' gout = annotate_C_emb(gout, "celltype")
#' names(colData(gout$C_emb))
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

#' propagate a colData variable to cell embedding; endomorphism
#' @param sscr instance of simba_scrna
#' @param orig_var character(1) must be element of `colData(sscr$ppCG)`
#' @return updated simba_scrna instance with `C_emb` colData updated
#' @export
annotate_C_emb = function(sscr, orig_var) {
 stopifnot(inherits(sscr, "simba_scrna"))
 orig = sscr$ppCG
 Cemb = sscr$C_emb
 vals = colData(orig)[[orig_var]]
 names(vals) = colnames(orig)
 vals_reord = vals[colnames(Cemb)]
 colData(Cemb)[[orig_var]] = as.character(vals_reord)
 sscr$C_emb = Cemb
 sscr
}

#' get a fitted simba_scrna for pbmc 3k
#' @export
get_fitted_3k = function() {
 readRDS(system.file("gout/gout3k.rds", package="BiocSIMBA"))
}

#' visualize the cell embedding
#' @import ggplot2
#' @importFrom uwot umap
#' @param gout instance of `simba_scrna` with annotated C_emb component
#' @param colour_by character(1) factor in colData of C_emb for colouring
#' @param ptsize for `geom_point` size
#' @param ltxtsize for `legend.text`
#' @param lkeysize for `legend.key.size` in cm
#' @param \dots passed to uwot::umap
#' @return a ggplot instance
#' @examples
#' g3k = get_fitted_3k()
#' viz_cemb_umap(g3k, "celltype")
#' @export
viz_cemb_umap = function(gout, colour_by, ptsize=5, ltxtsize=30,
   lkeysize=1.5, ...) {
  stopifnot(inherits(gout, "simba_scrna"))
  stopifnot(colour_by %in% names(colData(gout$C_emb)))
  um = uwot::umap(t(assay(gout$C_emb)), ...)
  dd = data.frame(x=um[,1], y=um[,2], fac=colData(gout$C_emb)[[colour_by]])
  ggplot(dd, aes(x=x,y=y,colour=fac, text=fac)) + geom_point(size=ptsize) + theme(legend.text=element_text(size=ltxtsize)) + theme(legend.key.size=unit(lkeysize, "cm"))
}


#ggplot(nn, aes(x=x,y=y,colour=type)) + geom_point(size=5) + theme(legend.text=element_text(size=30)) + theme(legend.key.size=unit(1.5, "cm"))

#' mutual embedding
#' @examples
#' g3k = get_fitted_3k()
#' jj = joint_emb_CG( g3k )
#' viz_joint_umap( jj, "celltype" )
#' @export
joint_emb_CG = function( sscr )
{
    proc = basilisk::basiliskStart(bsklenv, testload = "simba")
    basilisk::basiliskRun(proc, function( sscr ) {
     adC = zellkonverter::SCE2AnnData( sscr$C_emb )
     adG = zellkonverter::SCE2AnnData( sscr$G_emb )
     sref = reticulate::import("simba", convert=FALSE)
     adall = sref$tl$embed( adata_ref = adC, list_adata_query = list(adG) )
     zellkonverter::AnnData2SCE(adall)
    }, sscr )
}

#' view mutual embedding
#' @param scej a SingleCellExperiment instance produced by `joint_emb_CG`
#' @param cname character(1) an element of colData to be used for coloring
#' @param ptsize for `geom_point` size
#' @param ltxtsize for `legend.text`
#' @param lkeysize for `legend.key.size` in cm
#' @param \dots passed to uwot::umap
#' @return a ggplot instance
#' @export
viz_joint_umap = function (scej, cname, ptsize = 3, 
    ltxtsize = 30, lkeysize = 1.5, ...) {
    stopifnot(inherits(scej, "SingleCellExperiment"))
    stopifnot(cname %in% names(colData(scej)))
    uu = uwot::umap(t(assay(scej)), ...)
    txt = lab = colData(scej)[[cname]]
    if (any(isg <- is.na(lab))) 
        lab[is.na(lab)] = "gene"
    txt[which(isg)] = colnames(scej)[which(isg)]
    ndf = data.frame(x = uu[, 1], y = uu[, 2], fac = lab)
    ggplot(ndf, aes(x = x, y = y, colour = fac, text = txt)) + 
        geom_point(size = ptsize, alpha=.8) + theme(legend.text = element_text(size = ltxtsize)) + 
        theme(legend.key.size = unit(lkeysize, "cm"))
}
