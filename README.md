## R package of utility functions for ggplotting bio3d matrix data. 

The **bio3d.ggplot** package currently provides two categories of functions:
`gg.*()` and `gg_*()` functions. This first produce ggplot objects from
bio3d package related data.  The later add layers to the former.  

Currently the main function is `gg.mat()`, which returns a
ggplot version of the old bio3d `plot.dmat()` function for NxN matrix data. 
This function also works well for difference distance matrices (DDMs), 
correlation matrices (substituting for `plot.dccm()`) and contact maps 
(substituting for `plot.cmap()`).

Also included is a basic secondary structure element annotation function
called `gg_sse()` that provides the classic alpha helix and beta
strand rectangles seen in bio3d plot functions since day one. Options for 
more fancy helices will follow.  

This a work in progress.
As things develop I will likely introduce S3 based classes that call this one 
master matrix plot function with suitable settings. As well as one master 
vector plot function (the later like `plot.bio3d()`. The the `gg_*()` 
annotation functions should work for all of these in the same way.   

To install:

```
devtools::install_github("bioboot/bio3d.ggplot")
```

To make some nice figures (i.e. nicer than base graphics versions):

```
library(bio3d.ggplot)
library(bio3d)
library(ggplot2)

## Run the examples below
example(bio3d.ggplot)
```

Or run each example in turn (these are in the man/help page for the package):

```
 ##- Single structure distance matrix
 pdb <- read.pdb( "5p21" )
 k <- dm(pdb, inds="calpha", mask.lower=FALSE)
 gg.mat(k)

 ## Add secondary structure 'annotation layer' to plot
 p <- gg.mat(k) + gg_sse(pdb)  ## save in an object 'p'
 p  ## produce the plot


 ## Data driven axis from SSE boundaries
 sse_labels <- c(pdb$helix$start, pdb$sheet$start,
                 pdb$helix$end, pdb$sheet$end)
 x <- sort(sse_labels); x <- x[diff(x)> 4]

 ## Add your customization with additional layers
 p + theme_grey() + ## The ggplot2 gray background theme
  scale_fill_gradient(high = "orange", low = "white") +
  scale_y_continuous("My Y label",breaks=x, labels=x, expand=c(0,0.5)) +
  scale_x_continuous("My X label",breaks=x, labels=x, expand=c(0,0.5))


 gg.mat(k) + theme_minimal(base_size = 14) +
   scale_x_continuous(expand = c(0, 0)) +
   scale_y_continuous(expand = c(0, 0))  ## start axis at zero

 ##- Calculate and plot correlation matrix
 cij <- dccm.nma(nma(pdb))
 q <- gg.mat(cij) + gg_sse(pdb)
 q

 q +  scale_fill_gradient2(limit = c(-1,1),
          high = "red", mid = "white", low = "blue")

 ##- Difference distance matrices (DDM) of heterogeneous structures
 pdbs <- pdbaln( c("5p21","4q21") )
 mat <- dm(pdbs, mask.lower=FALSE)
 ddm <- mat[,,1] - mat[,,2]
 gg.mat(ddm) + gg_sse(pdbs)


 ##- Contact maps
 cm <- cmap(pdb, scut=0, mask.lower=FALSE)
 gg.mat(cm) + gg_sse(pdb) + theme(legend.position="none")
```
