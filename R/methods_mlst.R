#' @export
print.mlst <- function(x, ...){

  infiles <- attr(x, 'infiles')
  org <- attr(x, 'org')
  scheme <- attr(x, 'scheme')
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
#' @description Plot a \code{mlst} object. A Minimum Spanning Tree or a binary
#' tree is rendered.
#' @param x An object of class \code{mlst}.
#' @param type One of "mst" or "phylo", for plotting a Minimum Spanning Tree or
#' a binary tree.
#' @param what One of "result", "profile", "both" (default). What should be
#' plotted.
#' @param pt.size The size of the point. Default: 3.
#' @param label \code{logical}. Whether to plot node/tip labels or not. Default:
#' \code{FALSE}.
#' @param pf.col The color of profile nodes/tips. Ignored if \code{what="result"}.
#' @param st.col The color of result nodes/tips which have an ST assigned.
#' Ignored if \code{what="profile"}.
#' @param nst.col The color of result nodes/tips with no ST assigned. Default:
#' \code{"white"}. Ignored if \code{what="profile"}.
#' @param alpha Color transparency, between [0,1]. See \link[scales]{alpha}.
#' @param pairwise.deletion A \code{logical} indicating whether to delete the
#' comumns with missing data on a pairwise basis. If \code{FALSE} (default) and
#' \code{NA}s are found, then each \code{NA} is condidered as a different
#' allele when computing distance.
#' @param plot.igraph.args A \code{list} of arguments to be passed to
#' \link[igraph]{plot.igraph}. Used only if \code{type="mst"}. Defaults try
#' to keep aesthetics similar to \code{phylo}'s. Default: \code{list(vertex.label = if (label) NULL else NA,
#' vertex.size = pt.size)}.
#' @param plot.phylo.args  A \code{list} of arguments to be passed to
#' \link[ape]{plot.phylo}. Used only if \code{type="phylo"}. Defaults try to
#' keep aesthetics similar to \code{mst}'s. Default: \code{list(type='unrooted',
#' show.tip.label = label)}.
#' @param tiplabels.args  A \code{list} of arguments to be passed to
#' \link[ape]{tiplabels}. Used only if \code{type="phylo"}. Defaults try to
#' keep aesthetics similar to \code{mst}'s. Default: \code{list(pch = 21,
#' cex = pt.size/5)}.
#' @param plot Default: TRUE.
#' @param ... A list of arguments to be passed to \link[graphics]{par}.
#' @return A minimum spanning tree plot and an object of class \code{igraph}
#' (invisible) if \code{type="mst"}, or a binary tree plot and an object of
#' class \code{phylo} (invisible) if \code{type="phylo"}.
#' @details Distance is calculated using \link[ape]{dist.gene} function over
#' the allele matrix. This distance metric only takes into account the number
#' of differences between each pair of rows in the matrix. If \code{type="mst"},
#' \link[ape]{mst} function is used to calculate a minimum spanning tree, and
#' a graph is generated using \link[igraph]{graph.adjacency}. If
#' \code{type="phylo"}, a "phylogeny" is inferred using the distance calculated
#' above and \link[ape]{nj} function (neighbour-joining tree). Ethier a
#' \code{igraph} or a \code{phylo} object is returned invisibly so it can be
#' further analysed using igraph/ape frameworks.
#'
#' It is worth noting that the result of this function is not strictly based on
#' genetic information. The distance used just counts the number of differences
#' on the allele profiles, it doesn't use a genetic model. The clustering
#' methods are also simple and are not based on complex evolution models.
#' Despite the above, since alleles in mlst schemes usually differ in 1 or few
#' nucleotides, the described methodology gives good enough results to have a
#' general overview of your data.
#' @importFrom ape dist.gene mst nj plot.phylo tiplabels
#' @importFrom igraph graph.adjacency V<- V plot.igraph
#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics plot par
#' @importFrom scales alpha
#' @export
plot.mlst <- function(x,
                      type = "mst",
                      what = "both",
                      pt.size = 3,
                      label = FALSE,
                      pf.col = "#E64B35FF",
                      st.col = "#00A087FF",
                      nst.col = "white",
                      alpha = 0.5,
                      pairwise.deletion = FALSE,
                      plot.igraph.args = list(),
                      plot.phylo.args = list(),
                      tiplabels.args = list(),
                      plot = TRUE,
                      ...){

  type <- match.arg(type, choices = c('mst', 'phylo'))
  what <- match.arg(what, choices = c('result', 'profile', 'both'))

  if(class(label)!='logical'){
    stop('label must be logical.')
  }
  if(length(label)!=1L){
    stop('label must be a single value of class "logical".')
  }

  if(what=='result'){

    resu <- x$result
    nas <- is.na(resu$ST)
    sts <- resu$ST[which(!nas)]
    nst <- rownames(which(nas))
    di <- dim(resu)
    m <- resu[, -di[2]]
    cols <- rep(pf.col, dim(m)[1])
    cols[nst] <- nst.col

  }else if(what=='profile'){

    prof <- x$profile
    di2 <- dim(prof)
    m <- prof[, -di2[2]]
    cols <- rep(pf.col, dim(m)[1])

  }else{

    resu <- x$result
    nas <- is.na(resu$ST)
    sts <- which(!nas)
    nst <- which(nas)
    di <- dim(resu)
    resu2 <- resu[, -di[2]]

    prof <- x$profile
    di2 <- dim(prof)
    prof2 <- prof[, -di2[2]]

    m <- rbind(resu2, prof2)
    cols <- rep(pf.col, dim(m)[1])
    cols[sts] <- st.col
    cols[nst] <- nst.col
  }


  cols <- alpha(cols, alpha = alpha)

  # d <- dist.gene(m)/(dim(m)[2]-1L)
  # following is better cause it return integer and not numeric, and in
  # practical cases it gives the same result on further steps:
  nam <- is.na(m)
  if (any(nam) & !pairwise.deletion) {
    warning("NAs encountered.\n As pairwise.deletion = FALSE, each NA will be considered as a different allele.", immediate. = TRUE, noBreaks. = TRUE)
    nna <- length(which(nam))
    m[which(nam, arr.ind = TRUE)] <- paste(m[which(nam, arr.ind = TRUE)], seq_len(nna), sep = "")
  }
  d <- dist.gene(m, pairwise.deletion = pairwise.deletion)


  ### Start device functions ###
  on.exit(dev.flush())
  dev.hold()

  apar <- list(mar=c(0,0,0,0)+.1)
  inargs <- list(...)
  apar[names(inargs)] <- inargs
  do.call('par', apar)

  if (type=='mst'){

    tree <- mst(d)
    g <- graph.adjacency(tree, mode = 'undirected')
    V(g)$color <- cols

    if(plot){

      piargs <- list(vertex.label = if (label) NULL else NA,
                     vertex.size = pt.size)

      piargs[names(plot.igraph.args)] <- plot.igraph.args

      args <- c(list(g), piargs)

      do.call('plot', args)

    }

  }else{

    g <- nj(d)

    if(plot){

      ppargs <- list(type='unrooted',show.tip.label = label)
      ppargs[names(plot.phylo.args)] <- plot.phylo.args
      args1 <- c(list(g), ppargs)

      tlargs <- list(pch = 21,cex = pt.size/4)
      tlargs[names(tiplabels.args)] <- tiplabels.args
      args2 <- c(list(bg = cols), tlargs)

      do.call('plot', args1)
      do.call('tiplabels', args2)

    }

  }


  invisible(g)

}
