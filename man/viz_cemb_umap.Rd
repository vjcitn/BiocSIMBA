% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_graph_sce.R
\name{viz_cemb_umap}
\alias{viz_cemb_umap}
\title{visualize the cell embedding}
\usage{
viz_cemb_umap(gout, colour_by, ptsize = 5, ltxtsize = 30, lkeysize = 1.5, ...)
}
\arguments{
\item{gout}{instance of `simba_scrna` with annotated C_emb component}

\item{colour_by}{character(1) factor in colData of C_emb for colouring}

\item{ptsize}{for `geom_point` size}

\item{ltxtsize}{for `legend.text`}

\item{lkeysize}{for `legend.key.size` in cm}

\item{\dots}{passed to uwot::umap}
}
\value{
a ggplot instance
}
\description{
visualize the cell embedding
}
\examples{
g3k = get_fitted_3k()
viz_cemb_umap(g3k, "celltype")
}
