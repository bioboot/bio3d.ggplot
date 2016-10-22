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
#' @param side A numeric vector with values 1 to 4 specifying he side of the 
#'   plot where annotations will be drawn. Note sides 3 and 4 are not set 
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
#'  pdb <- read.pdb( "5p21" )
#'  k <- dm(pdb, inds="calpha", mask.lower=FALSE)
#'  gg.mat(k)
#'
#'  gg.mat(k) + gg_sse(pdb) ## Add secondary structure from pdb
#'
#' @export
gg_sse <- function(x, min=-5, max=0, helix.col="gray20", sheet.col="gray80", sse.border="black", side=c(1,2)) {
  ## Add secondary structure to a ggplot
  ##
  ## input x can be a 'pdb' class, 'pdbs' class or vector as obtained from pdb2sse()
  ## - ToDo: Add a side=1:4 option (currently always 1 and 2)
  ##         Remove redundancy in sse annotation calls for diff inputs below
  ##         Find absolute coords/scale for 'min' independent of data/plot coords...
  ##           Better yet have a 'sse.scale' and optional 'sse.size' arguments to control width.


  ### min=-5; max=0; helix.col="gray20"; sheet.col="gray80"; sse.border="black"; side=3
  #ymax=100+abs(min)
  #xmax=100+abs(min)
  ### Need to work out max plot dims in a better way!!
  if(!(is.numeric(side) && all(side %in% 1:4))) {
    stop("The 'side=' for SSE annotation should be a vector of values between 1 and 4 only")
  }

  mid <- (min-max)/2

  ##- Single structure input all sides will be the same.
  if( inherits(x, c("pdb","sse")) ) {

    if( inherits(x,"pdb") ) {

      ## Maximum plot dimension 'xmax' and 'ymax'
      ## For now we get max plot dim from residue number.
      ## NEED A BETTER WAY (i.e. either pass plot object or inheret it)
      xmax <- sum(summary(pdb)[c("nprot.res", "nother.res")])

      ## PDB objects dont currently rtn $sse vector from which to get nres...
      xend   <- sum(x$calpha)
    } else {
      xend   <- length(x$sse)
      xmax <- xend + 4 ## !!!-- Hack: we dont have non-protein in dssp output --!!!
    }

    ## For sides 1 and 3
    hstart <- x$helix$start
    hend   <- x$helix$end
    sstart <- x$sheet$start
    send   <- x$sheet$end

    ## For sides 2 and 4 
    ##  (one input structure only so these will be the same)
    hstart.2 <- hstart
    hend.2   <- hend
    sstart.2 <- sstart
    send.2   <- send

  }


  ## Two structure input assumed e.g. DDM
  if( inherits(x,"pdbs") ){

    xend <- ncol(x$sse)
    xmax <- xend

    h <- which(x$sse == "H", arr.ind=TRUE)
    e <- which(x$sse == "E", arr.ind=TRUE)

    ##-- N.B.-- Consider first two structures only!!
    h1 = bio3d::bounds(h[h[,1]==1,2], pre.sort = FALSE)
    e1 = bio3d::bounds(e[e[,1]==1,2], pre.sort = FALSE)

    h2 = bio3d::bounds(h[h[,1]==2,2], pre.sort = FALSE)
    e2 = bio3d::bounds(e[e[,1]==2,2], pre.sort = FALSE)

    ## For sides 1 and 3
    hstart <- h1[,"start"]
    hend   <- h1[,"end"]
    sstart <- e1[,"start"]
    send   <- e1[,"end"]

    ## For sides 2 and 4 
    ##  (note different structure here!)
    hstart.2 <- h2[,"start"]
    hend.2   <- h2[,"end"]
    sstart.2 <- e2[,"start"]
    send.2   <- e2[,"end"]


  } else {

    if(is.vector(x)) {
      ## Take a vector of input as produced from pdb2sse(pdb)
      h <- bio3d::bounds(which(x == "H"), pre.sort = FALSE)
      e <- bio3d::bounds(which(x == "E"), pre.sort = FALSE)

      xend <- length(x)
      xmax <- xend
      
      hstart <- h[,"start"]
      hend   <- h[,"end"]
      sstart <- e[,"start"]
      send   <- e[,"end"]
      hstart.2 <- hstart
      hend.2   <- hend
      sstart.2 <- sstart
      send.2   <- send
    }
  }

  ## Build up vectors for helix and sheet rectangles
  hxmin=NULL; hxmax=NULL; hymin=NULL; hymax=NULL
  sxmin=NULL; sxmax=NULL; symin=NULL; symax=NULL
  segx=NULL; segxend=NULL; segy=NULL; segyend=NULL

  # Plot max limit
  xmax=xmax+abs(min)
  ymax=xmax

  # side1
  if( any(side == 1) ) {
    segx <- c(segx, 1); segxend <- c(segxend, xend)
    segy <- c(segy, mid); segyend <- c(segyend, mid)

    hxmin <- c(hxmin, hstart) 
    hxmax <- c(hxmax, hend)
    hymin <- c(hymin, rep(min, length(hstart)))
    hymax <- c(hymax, rep(max, length(hstart)))

    sxmin <- c(sxmin, sstart) 
    sxmax <- c(sxmax, send)
    symin <- c(symin, rep(min, length(sstart)))
    symax <- c(symax, rep(max, length(sstart)))
  }

  # side2
  if( any(side == 2) ) {
    segx <- c(segx, mid); segxend <- c(segxend, mid)
    segy <- c(segy, 1); segyend <- c(segyend, xend)

    hxmin <- c(hxmin, rep(min, length(hstart))) 
    hxmax <- c(hxmax, rep(max, length(hstart)))
    hymin <- c(hymin, hstart.2)
    hymax <- c(hymax, hend.2)

    sxmin <- c(sxmin, rep(min, length(sstart))) 
    sxmax <- c(sxmax, rep(max, length(sstart)))
    symin <- c(symin, sstart.2)
    symax <- c(symax, send.2)
  }

  # side3
  if( any(side == 3) ) {
    segx <- c(segx, 1); segxend <- c(segxend, xend)
    segy <- c(segy, mid+ymax); segyend <- c(segyend, mid+ymax)

    hxmin <- c(hxmin, hstart) 
    hxmax <- c(hxmax, hend)
    hymin <- c(hymin, rep((min+ymax), length(hstart)))
    hymax <- c(hymax, rep((max+ymax), length(hstart)))

    sxmin <- c(sxmin, sstart) 
    sxmax <- c(sxmax, send)
    symin <- c(symin, rep((min+ymax), length(sstart)))
    symax <- c(symax, rep((max+ymax), length(sstart)))
  }

  # side4
  if( any(side == 4) ) {
    segx <- c(segx, mid+xmax); segxend <- c(segxend, mid+xmax)
    segy <- c(segy, 1); segyend <- c(segyend, xend)

    hxmin <- c(hxmin, rep((min+xmax), length(hstart))) 
    hxmax <- c(hxmax, rep((max+xmax), length(hstart)))
    hymin <- c(hymin, hstart.2)
    hymax <- c(hymax, hend.2)

    sxmin <- c(sxmin, rep((min+xmax), length(sstart))) 
    sxmax <- c(sxmax, rep((max+xmax), length(sstart)))
    symin <- c(symin, sstart.2)
    symax <- c(symax, send.2)
  }

  out <- list(
    ggplot2::annotate("segment", x=segx, xend=segxend, y=segy, yend=segyend),
    ggplot2::annotate("rect", xmin=hxmin, xmax=hxmax, ymin=hymin, ymax=hymax, col=sse.border, bg=helix.col),
    ggplot2::annotate("rect", xmin=sxmin, xmax=sxmax, ymin=symin, ymax=symax, col=sse.border, bg=sheet.col)
  )
}

