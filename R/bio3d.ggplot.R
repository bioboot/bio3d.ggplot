#' bio3d.ggplot: Easy 'ggplotting' of bio3d package Data.
#'
#' The bio3d.ggplot package provides two categories of functions:
#'  ggplot.* and gg_* functions. This first produce the ggplot objects from 
#'  bio3d data.  The later add layers to the former. 
#'
#' @note Currently the main function is \code{ggplot.dmat}, which returns a 
#'  ggplot version of the old plot.dmat() for NxN matrix data.
#'
#'  Also included is a basic secondary structure element annotation function 
#'  called \code{gg_sse()} that provides the classic alpha helix and beta 
#'  strand rectangles seen in bio3d since day one. Options for more fancy  
#'  helices will follow.
#'
#'  This a work in progress.
#'  As things develop I will likely have one master matrix plot function here
#'  and one vector plot function (the later like \code{plot.bio3d}. The 
#'  plot.dmat function is the start of the first. 
#'
#' @docType package
#'
#' @name bio3d.ggplot
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
#'
#' \donttest{ 
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
#'
#'  ##- Calculate and plot correlation matrix
#'  cij <- dccm.nma(nma(pdb))
#'  q <- ggplot.dmat(cij) + gg_sse(pdb)
#'  q
#'
#'  q +  scale_fill_gradient2(limit = c(-1,1),
#'           high = "red", mid = "white", low = "blue")
#'
#'
#'  ##- Difference distance matrices (DDM) of heterogeneous structures
#'  pdbs <- pdbaln( c("5p21","4q21") )
#'  mat <- dm(pdbs, mask.lower=FALSE)
#'  ddm <- mat[,,1] - mat[,,2]
#'  ggplot.dmat(ddm) + gg_sse(pdbs)
#'
#'
#'  ##- Contact maps
#'  cm <- cmap(pdb, scut=0, mask.lower=FALSE)
#'  ggplot.dmat(cm) + gg_sse(pdb) + theme(legend.position="none")
#' }
#' 
NULL
