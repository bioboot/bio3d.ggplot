#' ggplot bio3d matrix data
#'
#' \code{ggplot.dmat} returns a ggplot version of plot.dmat() NxN matrix plot.
#'
#' This is a simple first go at getting our usual bio3d matrices into a ggplot
#'   friendly format. This works with utility functions for adding annotations
#'   such as secondary structure, domain boundarys highlight regions etc.
#'
#' @param dm A numeric matrix such as a distance matrix from the bio3d functions
#' \code{dm}, \code{cmap}, \code{dccm} etc. Could write these as \code{\link[bio3d]{dm}},
#' \code{\link[bio3d]{cmap}}, \code{\link[bio3d]{dccm}}, \code{pca.array.loadings},
#' etc.
#'
#' @return A ggplot object for plotting with \code{print.ggplot}.
#'
#' @note Still much to be improved here, including adding extra args via dots.
#'
#' @seealso The bio3d base, grid and lattice ploting functions for various matrix class objects.
#'   This includes \code{\link[bio3d]{plot.dmat}}, \code{\link[bio3d]{plot.dccm}},
#'   \code{\link[bio3d]{plot.cmap}}, etc.
#'
#' @examples
#'  pdb <- bio3d::read.pdb( "5p21" )
#'  k <- bio3d::dm(pdb, inds="calpha", mask.lower=FALSE)
#'  ggplot.dmat(k)
#'
#'  ggplot.dmat(k) + gg_sse(pdb) ## Add secondary structure from pdb
#'
#'  ggplot(k) + theme_grey() +   ## The ggplot2 gray background theme
#'    scale_x_continuous(expand = c(0, 0)) +
#'    scale_y_continuous(expand = c(0, 0))  ## start axis at zero

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

  p <- ggplot(data = h, aes(x=Var1, y=Var2, fill=value)) + geom_raster() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name="") +
    theme_bw(base_size = 14) +
    scale_x_continuous(expand = c(0, 0.1)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    labs(x="Residue Number", y="Residue Number") +
    coord_fixed()

  return(p)
}
