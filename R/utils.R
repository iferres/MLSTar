processResu <- function(resu,
                        write,
                        dbs,
                        dir,
                        dnw = paste0(dir, 'tmp/'),
                        prefix){

  #Cat files
  sapply(basename(dbs), function(x){
    ge <- list.files(path = dnw,
                     pattern = x,
                     recursive = TRUE,
                     full.names = TRUE)
    if (length(ge)>0){
      file.copy(ge[1], to = dnw)
      nw <- paste0(dnw, basename(ge[1]))
      if(length(ge)>1){
        for (i in 2:length(ge)){
          file.append(nw, ge[i])
        }
      }
    }
  })
  unlink(list.dirs(dnw)[-1], recursive = TRUE)

  if (write%in%c('new', 'all')){
    sapply(basename(dbs), function(x){
      ge <- list.files(path = paste0(dir, prefix),
                       pattern = x,
                       recursive = TRUE,
                       full.names = TRUE)
      if (length(ge)>0){
        file.copy(ge[1], to = paste0(dir, prefix))
        nw <- paste0(paste0(dir, prefix), '/', basename(ge[1]))
        if(length(ge)>1){
          for (i in 2:length(ge)){
            file.append(nw, ge[i])
          }
        }
      }
    })
    unlink(list.dirs(paste0(dir, prefix))[-1], recursive = TRUE)
  }


  #Check if 'u'nknowns are not the same, and name accordingly.
  us <- which(apply(resu, 2,function(x){'u'%in%x}))
  count <- 1L
  if (length(us)>0){
    cat('Checking if new alleles are equal...')
    for (u in seq_along(us)){
      gene <- names(us[u])
      out.tmp <- paste0(dnw, gene, '.fasta')
      gp <- grep('u', resu[, us[u]])
      if (length(gp)>1){
        sq <- seqinr::read.fasta(out.tmp, seqtype = 'DNA', as.string = TRUE)
        ul <- unlist(sq)
        gr <- split(names(ul), ul)
        for (i in seq_along(gr)){
          nu <- paste0('u', count)
          n <- sapply(strsplit(gr[[i]],';'), '[', 2)
          resu[n, us[u]] <- nu
          if (write%in%c('new','all')){
            out.newAllele <- paste0(dir, prefix, '/', gene, '.fasta')
            sq <- seqinr::read.fasta(out.newAllele, seqtype = 'DNA', as.string = TRUE)
            nsq <- sapply(paste0(';', n, ';'), function(z){ grep(z, names(sq), fixed = TRUE)})
            spl <- do.call(rbind, strsplit(names(sq)[nsq], ';'))
            spl[, 1] <- sub('u$', nu, spl[, 1])
            names(sq)[nsq] <- paste0(spl, collapse = ';')
            seqinr::write.fasta(sq,names = names(sq), file.out = out.newAllele, as.string = TRUE)
          }
          count <- count +1L
        }
      }else{
        nu <- paste0('u', count)
        resu[gp, us[u]] <- nu
        if (write%in%c('new','all')){
          out.newAllele <- paste0(dir, prefix, '/', gene, '.fasta')
          is <- sub('[.]\\w+$','',rownames(resu)[gp])
          sq <- seqinr::read.fasta(out.newAllele, seqtype = 'DNA', as.string = TRUE)
          nsq <- grep(is, names(sq), fixed = TRUE)
          spl <- strsplit(names(sq)[nsq], ';')[[1]]
          spl[1] <- sub('u$', nu, spl[1])
          names(sq)[nsq] <- paste0(spl, collapse = ';')
          seqinr::write.fasta(sq,names = names(sq), file.out = out.newAllele, as.string = TRUE)
        }
        count <- count +1L
      }
    }
    cat(' DONE!\n')
  }

  unlink(dnw, recursive = TRUE)

  resu
}


#' @name mlst
#' @title MLST
#' @description Perform Multi Locus Sequence Typing.
#' @param genome A FASTA file with the genome sequence.
#' @param dbs Subject database files. No subfix needed.
#' @param write \code{character}. One of \code{"new"} (Default), \code{"all"} or
#' \code{"none"}. The fist one writes only new alleles found (not reported in
#' pubmlst.org), the second writes all alleles found, and "none" do not write
#' any file.
#' @param prefix \code{character} A prefix to the fasta files of found
#' sequences (see \code{write}). (Default: \code{"allele"}).
#' @param dir An existing directory where to put the loci fasta files in case
#' they are not provided by the user. Also sequences found will be placed here
#' if \code{write} is set to ethier \code{"new"} or \code{"all"}.
#' @param dnw Temporary directory where to put new alleles to check later if
#' are the same.
#' @param n_threads \code{integer}. The number of threads to use by BLASTN.
#' @param outf Where the blastn output will be written.
#' @param pid Percentage identity.
#' @param scov Query coverage.
#' @return A vector with the mlst.
#' @author Ignacio Ferres
#' @references  Altschul, Gish, Miller, Myers & Lipman. (1990) "Basic local
#' alignment search tool." \emph{J. Mol. Biol}. \strong{215}:403-410.
#'
#' Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming interface
#' for the PubMLST molecular typing and genome databases. Database: The Journal
#' of Biological Databases and Curation. 2017;2017:bax060.
#' doi:10.1093/database/bax060.
#' @importFrom parallel mclapply
#' @importFrom utils read.csv
mlst <- function(genome,
                 dbs,
                 write='new',
                 prefix = 'allele',
                 dir='.',
                 dnw = paste0(dir, 'tmp/'),
                 n_threads=1L,
                 outf=tempdir(),
                 pid=90,
                 scov=0.9){

  paste0(normalizePath(dir),'/') -> dir

  res <- c()
  nams <- c()
  for (i in 1:length(dbs)){
    # outf <- tempdir()
    blout <- blastn(genome = genome,
                    db = dbs[i],
                    outdir = outf,
                    n_threads = n_threads)
    blastRes <- readBlastResult(blout = blout)
    a <- processBlastResult(blastRes = blastRes,
                            write=write,
                            prefix = prefix,
                            dir = dir,
                            dnw = dnw,
                            pid = pid,
                            scov = scov)
    res[i] <- as.character(a)
    nams[i] <- names(a)
  }

  names(res) <- nams
  res
}


