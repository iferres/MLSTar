#' @export
print.mlst <- function(x, ...){

  infiles <- attr(mlst, 'infiles')
  org <- attr(mlst, 'org')
  scheme <- attr(mlst, 'scheme')
  pr <- paste0('MLST\n ',
               length(infiles),
               ' isolates\n organism: ',
               org,
               '\n scheme: ',
               scheme,
               '\n')
  cat(pr)

}


# as.data.frame.mlst <- function(mlst){
#   mlst$result
# }

# summary.mlst <- function(mlst){
#
#
#
# }

#' @name plot.mlst
#' @title Plot A mlst Object
#' @description Plot a \code{mlst} object. A Minimum Spanning Tree is generated
#' and a graph plot is rendered.
#' @param x An object of class \code{mlst}.
#' @param what One of "result", "profile", "both" (default). What should be
#' plotted. White nodes are plotted for those isolates with no ST assigned. If
#' "both", isolates with assigned ST are plotted in blue.
#' @param vertex.size The size of the vertex. Default: 3.
#' @param vertex.label Default: NA.
#' @param plot Default: TRUE.
#' @param ... Further arguments to pass to \link[igraph]{plot.igraph}.
#' @return A minimum spanning tree plot and an object of class \code{igraph}
#' (invisible).
#' @importFrom ape dist.gene mst
#' @importFrom igraph graph.adjacency V<- V plot.igraph
#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics plot par
#' @export
plot.mlst <- function(x,
                      what = 'both',
                      vertex.size = 3,
                      vertex.label = NA,
                      plot = TRUE,
                      ...){

  what <- match.arg(what, choices = c('result', 'profile', 'both'))

  if(what=='result'){

    resu <- x$result
    nas <- is.na(resu$ST)
    sts <- resu$ST[which(!nas)]
    nst <- rownames(which(nas))
    di <- dim(resu)
    m <- resu[, -di[2]]

  }else if(what=='profile'){

    prof <- x$profile
    di2 <- dim(prof)
    m <- prof[, -di2[2]]

  }else{

    resu <- x$result
    nas <- is.na(resu$ST)
    sts <- resu$ST[which(!nas)]
    nst <- rownames(which(nas))
    di <- dim(resu)
    resu2 <- resu[, -di[2]]

    prof <- x$profile
    di2 <- dim(prof)
    prof2 <- prof[, -di2[2]]

    m <- rbind(resu2, prof2)

  }


  d <- dist.gene(m)
  tree <- mst(d)
  g <- graph.adjacency(tree, mode = 'undirected')


  ### Start device functions ###
  op <- par(no.readonly = T)
  on.exit(dev.flush())
  on.exit(par(op), add = TRUE)
  dev.hold()


  if (what=='result'){
    V(g)$color <- 1
    V(g)[nst]$color <- NA
  }else if(what=='both'){
    V(g)$color <- 1
    V(g)[sts]$color <- 2
    V(g)[nst]$color <- NA
  }


  if(plot){

    par(mar=c(0,0,0,0)+.1)
    plot(g,
         vertex.label = vertex.label,
         vertex.size = vertex.size,
         ...)

  }


  invisible(g)

}
