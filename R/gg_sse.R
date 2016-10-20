#' Add Secondary Structure Annotation to a ggplot Matrix
#'
#' Add marginal secondary structure elements (SSEs) annotation layer to the
#'   current ggplot object.
#'
#' @param x A bio3d object with secondary structure information. Most commonly
#'   this is a pdb or pdbs class object as obtained from the bio3d functions
#'   read.pdb and read.fasta.pdbs. However, a vector of sse elements as obtained
#'   from the function pdb2sse is also valid input. The later is useful for more
#'   complicated multi-chain input structures. See examples.
#' @param min Currently the bottom most coordinate  for annotation rectangles. 
#'   More negative values will result in larger width rectangles. Eventually
#'   we will change this to introduce a single input for specifying the relative
#'   width of sse annotation rectangles.
#' @param max the top most plot coordinate for annotation rectangles.
#' @param helix.col The fill colors for rectangles representing alpha helices.
#' @param sheet.col The fill colors for rectangles representing beta strands.
#' @param sse.border The stroke color for rectangle borders.
#'
#' @return A list of annotation layers for adding to a ggplot object.
#'
#' @note There is lots to clean up and improve in this function. This includes 
#'   making on one call to annotate() for all rectangles rather than the four  
#'   calls per input type we do currently. Also we need a side=c(1:4) option  
#'   so we can add this for 1d plots and 2d plots (with the later providing any
#'    combination of annotated axis). 
#'
#'   Two forms of secondary structure annotation are available: so
#'   called ‘classic’ and ‘fancy’. The former draws marginal rectangles
#'   and has been available within Bio3D from version 0.1. The later
#'   draws more ‘fancy’ (and distracting) 3D like helices and arrowed strands.
#'   Also ToDo: add side options to allow drawing on sides 1 to 4
#'
#'
#' @examples
#'  pdb <- bio3d::read.pdb( "5p21" )
#'  k <- bio3d::dm(pdb, inds="calpha", mask.lower=FALSE)
#'  ggplot.dmat(k)
#'
#'  ggplot.dmat(k) + gg_sse(pdb) ## Add secondary structure from pdb
#'
#' @export

gg_sse <- function(x, min=-5, max=0, helix.col = "gray20", sheet.col = "gray80", sse.border = "black") {
  ## Add secondary structure to a ggplot
  ##
  ## input x can be a 'pdb' class, 'pdbs' class or vector as obtained from pdb2sse()
  ## - ToDo: Add a side=1:4 option (currently always 1 and 2)
  ##         Remove redundancy in sse annotation calls for diff inputs below
  ##         Find absolute coords/scale for 'min' independent of data/plot coords...
  ##           Better yet have a 'sse.scale' and optional 'sse.size' arguments to control width.


pdbs2helix <- function(pdbs) {
  ##- Convert 'sse' info from a 'pdbs' object to list format for gg_sse()
  h <- which(pdbs$sse == "H", arr.ind=TRUE)
  e <- which(pdbs$sse == "E", arr.ind=TRUE)

  ## Consider first two structures only!!
  h1 = bounds(h[h[,1]==1,2], pre.sort = FALSE)
  e1 = bounds(e[e[,1]==1,2], pre.sort = FALSE)

  h2 = bounds(h[h[,1]==2,2], pre.sort = FALSE)
  e2 = bounds(e[e[,1]==2,2], pre.sort = FALSE)

  ## could add a length filter here
  pdb1 <- list(helix=list(start=h1[,"start"], end=h1[,"end"]),
               sheet=list(start=e1[,"start"], end=e1[,"end"]) )

  pdb2 <- list(helix=list(start=h2[,"start"], end=h2[,"end"]),
               sheet=list(start=e2[,"start"], end=e2[,"end"]) )

  return( list(pdb1=pdb1, pdb2=pdb2 ) )
}



  mid <- (min-max)/2

  if( inherits(x,"pdb") ) {

    xend <- sum(x$calpha)
    out <- list(
      ggplot2::annotate("segment", x=1, xend=xend, y=mid, yend=mid),
      ggplot2::annotate("rect", xmin=x$helix$start, xmax =x$helix$end, ymin=min, ymax=max, col=sse.border, bg=helix.col),
      ggplot2::annotate("rect", xmin=x$sheet$start, xmax =x$sheet$end, ymin=min, ymax=max, col=sse.border, bg=sheet.col),

      ggplot2::annotate("segment", x=mid, xend=mid, y=1, yend=xend),
      ggplot2::annotate("rect", xmin=min, xmax=max, ymin=x$helix$start, ymax=x$helix$end, col=sse.border, bg=helix.col),
      ggplot2::annotate("rect", xmin=min, xmax=max, ymin=x$sheet$start, ymax=x$sheet$end, col=sse.border, bg=sheet.col) )
  }

  if( inherits(x,"pdbs") ){

    xs <- pdbs2helix(x)
    xend <- ncol(x$sse)
    out <- list(
      ggplot2::annotate("segment", x=1, xend=xend, y=mid, yend=mid),
      ggplot2::annotate("rect", xmin=xs$pdb1$helix$start, xmax=xs$pdb1$helix$end, ymin =min, ymax=max, col=sse.border, bg=helix.col),
      ggplot2::annotate("rect", xmin=xs$pdb1$sheet$start, xmax=xs$pdb1$sheet$end, ymin =min, ymax=max, col=sse.border, bg=sheet.col),

      ggplot2::annotate("segment", x=mid, xend=mid, y=1, yend=xend),
      ggplot2::annotate("rect", xmin=min, xmax=max, ymin=xs$pdb2$helix$start, ymax=xs$pdb2$helix$end, col=sse.border, bg=helix.col),
      ggplot2::annotate("rect", xmin=min, xmax=max, ymin=xs$pdb2$sheet$start, ymax=xs$pdb2$sheet$end, col=sse.border, bg=sheet.col) )

  } else {

    if(is.vector(x)) {
      ## Take a vector of input as produced from pdb2sse(pdb)
      h <- bounds(which(x == "H"), pre.sort = FALSE)
      e <- bounds(which(x == "E"), pre.sort = FALSE)

      xend <- length(x)
      out <- list(
        ggplot2::annotate("segment", x=1, xend=xend, y=mid, yend=mid),
        ggplot2::annotate("rect", xmin=h[,"start"], xmax=h[,"end"], ymin =min, ymax=max, col=sse.border, bg=helix.col),
        ggplot2::annotate("rect", xmin=e[,"start"], xmax=e[,"end"], ymin =min, ymax=max, col=sse.border, bg=sheet.col),

        ggplot2::annotate("segment", x=mid, xend=mid, y=1, yend=xend),
        ggplot2::annotate("rect", xmin=min, xmax=max, ymin=h[,"start"], ymax=h[,"end"], col=sse.border, bg=helix.col),
        ggplot2::annotate("rect", xmin=min, xmax=max, ymin=e[,"start"], ymax=e[,"end"], col=sse.border, bg=sheet.col) )
    }
  }
  out
}
