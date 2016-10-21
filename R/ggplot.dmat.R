#' ggplot bio3d matrix data
#'
#' \code{ggplot.dmat} returns a ggplot version of plot.dmat() NxN matrix plot.
#'
#' This is a simple first go at getting our usual bio3d matrices into a ggplot
#'   friendly format. This works with utility functions for adding annotations
#'   such as secondary structure, domain boundaries, highlight regions etc.
#'
#' @param dm A numeric matrix such as a distance matrix from the bio3d
#'   functions \code{dm}, \code{cmap}, \code{dccm} etc. Could write these
#'   as \code{\link[bio3d]{dm}}, \code{\link[bio3d]{cmap}}, \code{\link[bio3d]{dccm}},
#'   \code{pca.array.loadings}, etc.
#'
#' @return A ggplot object for plotting with \code{print.ggplot}.
#'
#' @note Still much to be improved here, including adding extra args via dots.
#'
#' @seealso The function \code{\link{gg_sse}} for adding secondary structure
#'   annotation. The various bio3d base, grid and lattice plotting functions
#'   for various bio3d matrix class objects. This includes:
#'  \code{\link[bio3d]{plot.dmat}},
#'   \code{\link[bio3d]{plot.dccm}}, \code{\link[bio3d]{plot.cmap}}, etc.
#'
#' @examples
#'  ##- Single structure distance matrix
#'  pdb <- bio3d::read.pdb( "5p21" )
#'  k <- bio3d::dm(pdb, inds="calpha", mask.lower=FALSE)
#'  ggplot.dmat(k)
#'
#'  ## Add secondary structure 'annotation layer' to plot
#'  p <- ggplot.dmat(k) + gg_sse(pdb)  ## save in an object 'p'
#'  p  ## produce the plot
#'
#' \dontrun{
#'  ## Data driven axis from SSE boundaries
#'  sse_labels <- c(pdb$helix$start, pdb$sheet$start,
#'                  pdb$helix$end, pdb$sheet$end)
#'  x <- sort(sse_labels); x <- x[diff(x)> 4]
#'
#'  ## Add your customization with additional layers
#'  p + theme_grey() + ## The ggplot2 gray background theme
#'   scale_fill_gradient(high = "orange", low = "white") +
#'   scale_y_continuous("My Y label",breaks=x, labels=x, expand=c(0,0.5)) +
#'   scale_x_continuous("My X label",breaks=x, labels=x, expand=c(0,0.5))
#'
#'
#'  ggplot.dmat(k) + theme_minimal(base_size = 14) +
#'    scale_x_continuous(expand = c(0, 0)) +
#'    scale_y_continuous(expand = c(0, 0))  ## start axis at zero
#'
#' ##- Calculate and plot correlation matrix
#' cij <- dccm.nma(nma(pdb))
#' q <- ggplot.dmat(cij) + gg_sse(pdb)
#' q
#'
#' q +  scale_fill_gradient2(limit = c(-1,1),
#'          high = "red", mid = "white", low = "blue")
#'
#' ##- Difference distance matrices (DDM) of heterogeneous structures
#' pdbs <- pdbaln( c("5p21","4q21") )
#' mat <- dm(pdbs, mask.lower=FALSE)
#' ddm <- mat[,,1] - mat[,,2]
#' ggplot.dmat(ddm) + gg_sse(pdbs)
#'
#' ##- Contact maps
#' cm <- cmap(pdb, scut=0, mask.lower=FALSE)
#' ggplot.dmat(cm) + gg_sse(pdb) + theme(legend.position="none")
#'
#' }
#' @export

ggplot.dmat <- function(dm){
  ## Simple ggplot version of plot.dmat()
  ##   ToDo: Add all the options as input args and dots() ...

  class(dm) = "matrix"

  if( is.null(rownames(dm)) ) {
    rownames(dm) = 1:nrow(dm)
  }

  if( is.null(colnames(dm)) ) {
    colnames(dm) = 1:ncol(dm)
  }

  h <- reshape2::melt(dm, na.rm = TRUE)

  p <- ggplot2::ggplot(data = h, ggplot2::aes_string(x="Var1", y="Var2", fill="value")) +
          ggplot2::geom_raster() +
          ggplot2::scale_fill_gradient2(low = "blue",
                                        high = "red",
                                        mid = "white",
                                        midpoint = 0,
                                        name="") +
          ggplot2::theme_bw(base_size = 14) +
          ggplot2::scale_x_continuous(expand = c(0, 0.1)) +
          ggplot2::scale_y_continuous(expand = c(0, 0.1)) +
          ggplot2::labs(x="Residue Number", y="Residue Number") +
          ggplot2::coord_fixed()

  return(p)
}
