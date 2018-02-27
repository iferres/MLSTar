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


as.data.frame.mlst <- function(mlst){
  mlst$result
}

summary.mlst <- function(mlst){



}


plot.mlst <- function(mlst,
                      vertex.size = 5,
                      vertex.label = NA,
                      layout = layout.fruchterman.reingold,
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

  d <- ape::dist.gene(m)
  tree <- ape::mst(d)
  g <- igraph::graph.adjacency(tree,
                               mode = 'undirected')

  V(g)$color <- 1
  V(g)[sts]$color <- 2
  V(g)[nst]$color <- NA

  plot(g,
       layout = layout,
       vertex.label = vertex.label,
       vertex.size = vertex.size,
       ...)

  invisible(g)

}
