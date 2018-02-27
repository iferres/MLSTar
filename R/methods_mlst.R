#' @export
print.mlst <- function(mlst){

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
#' @param mlst An object of class \code{mlst}.
#' @param vertex.size The size of the vertex. Default: 3.
#' @param vertex.label Default: NA.
#' @param layout layout. Default: layout.fruchterman.reingold.
#' @param plot Default: TRUE.
#' @param ... Further arguments to pass to \link[igraph]{plot.igraph}.
#' @return A minimum spanning tree plot and an object of class \code{igraph}
#' (invisible).
#' @importFrom ape dist.gene mst
#' @importFrom igraph graph.adjacency V plot.igraph
plot.mlst <- function(mlst,
                      vertex.size = 3,
                      vertex.label = NA,
                      layout = layout.fruchterman.reingold,
                      plot = TRUE,
                      ...){

  resu <- mlst$result
  nas <- is.na(resu$ST)
  sts <- resu$ST[which(!nas)]
  nst <- rownames(which(nas))
  di <- dim(resu)
  resu2 <- resu[, -di[2]]

  prof <- mlst$profile
  di2 <- dim(prof)
  prof2 <- prof[, -di2[2]]

  m <- rbind(resu2, prof2)

  d <- dist.gene(m)
  tree <- mst(d)
  g <- graph.adjacency(tree,
                               mode = 'undirected')


  ### Start device functions ###
  op <- par(no.readonly = T)
  on.exit(dev.flush())
  on.exit(par(op), add = TRUE)
  dev.hold()


  V(g)$color <- 1
  V(g)[sts]$color <- 2
  V(g)[nst]$color <- NA


  if(plot){

    par(mar=c(0,0,0,0)+.1)
    plot(g,
         layout = layout,
         vertex.label = vertex.label,
         vertex.size = vertex.size,
         ...)

  }


  invisible(g)

}
