#' @name doMLST
#' @title Perform MLST Analysis Over a List of Genomes
#' @description Takes a list of genome fasta files and perform blastn searches
#' to identify the sequence type for each of the genes/loci available in a mlst
#' scheme.
#' @param infiles A vector of genome sequences in fasta format.
#' @param org A valid organism from pubmlst.org. Run \link{listPubmlst_orgs} to
#' see available ones.
#' @param scheme \code{integer}. The scheme id number for a given organism. Run
#' \link{listPubmlst_schemes} to see available schemes for certain \code{org}.
#' @param schemeFastas A vector with the path to fasta sequences from each
#' loci of the specified mlst scheme. If it is \code{NULL} (default), sequences
#' are downloaded from \url{http://rest.pubmlst.org} to \code{dir}.
#' @param schemeProfile The path to the profile file (.tab). If left \code{NULL}
#' then it is downloaded from \url{http://rest.pubmlst.org} to \code{dir}.
#' @param write \code{character}. One of \code{"new"} (Default), \code{"all"} or
#' \code{"none"}. The fist one writes only new alleles found (not reported in
#' pubmlst.org), the second writes all alleles found, and "none" do not write
#' any file.
#' @param ddir A non-existing directory where to download the loci fasta files
#' in case they are not provided by the user. Default:
#' \code{paste0('pubmlst','_',org,'_',scheme)}
#' @param fdir A non-existing directory where to write fasta files of found
#' sequences (see \code{write}). (Default:
#' \code{paste0('alleles','_',org,'_',scheme)}).
#' @param n_threads \code{integer}. The number of cores to use. Each job consist
#' on the process for one genome. Blastn searches will use 1 core per job.
#' @param pid Percentage identity threshold to be consider as an allele. An
#' \code{integer} <= 100. (Default: 90).
#' @param scov Subject coverage threshold to be consider as an allele. A
#' \code{numeric} between 0 and 1. Not recomended to set it below 0.7 .
#' (Default 0.9)
#'# @details
#' @return A \code{data.frame}. Rows are \code{infiles} and columns are genes
#' from the selected MLST scheme. '<NA>' values means that no allele were
#' found in the respective genome. 'u' (from 'unknown') means that the allele
#' found was not reported in pubmlst.org database; a fasta file with the new
#' allele is written in this case. A number indicates the allele id number of
#' the reported alleles in pubmlst.org . The last column ($ST) indicates the ST
#' determined by the sequence of alleles, if available.
#' @author Ignacio Ferres
#' @references  Altschul, Gish, Miller, Myers & Lipman. (1990) "Basic local
#' alignment search tool." \emph{J. Mol. Biol}. \strong{215}:403-410.
#'
#' Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming interface
#' for the PubMLST molecular typing and genome databases. Database: The Journal
#' of Biological Databases and Curation. 2017;2017:bax060.
#' doi:10.1093/database/bax060.
#' @importFrom parallel mclapply
#' @export
doMLST <- function(infiles,
                   org='leptospira',
                   scheme=1L,
                   schemeFastas=NULL,
                   schemeProfile=NULL,
                   write = 'new',
                   ddir = paste0('pubmlst','_',org,'_',scheme),
                   fdir = paste0('alleles','_',org,'_',scheme),
                   n_threads=1L,
                   pid=90L,
                   scov=0.9){

  if (.Platform$OS.type!='unix'){
    stop('Sorry, this package works only on unix-like platforms.')
  }

  if(Sys.which('blastn')==''){
    stop('blastn binary is not in $PATH. Please install it before running this function.')
  }

  if (is.null(schemeFastas) | is.null(schemeProfile)){
    if(!dir.exists(ddir)){
      dir.create(ddir)
      ddir <- paste0(normalizePath(ddir),'/')
    }else{
      stop(paste0(ddir, ' already exists.'))
    }
  }

  write <- match.arg(write, c('none', 'new', 'all'))

  if (write%in%c('new', 'all')){
    if(!dir.exists(fdir)){
      #dont normalize path because in following steps directories are created
      # recursively.
      dir.create(fdir)
    }else{
      stop(paste0(fdir, ' already exists.'))
    }
  }

  if(pid>100){
    stop('pid must be a integer smaller than 100.')
  }

  if(scov>1 | scov<0){
    stop('scov must be between 0 and 1')
  }else if(scov < 0.7){
    warning('scov below 0.7 . Recomended: scov >= 0.7 . Continuing anyway..',
            immediate. = T)
  }

  if(any(!file.exists(infiles))){
    stop('One or more infiles do/es not exists in specified path.')
  }


  #Valid org?
  if(!org%in%listPubmlst_orgs()){
    stop(paste0(org,' is not in pubmlst database. Try listPubmlst_orgs() to see avaliables.'))
  }


  #Fasta files: If NULL, download; else, check if ok.
  if(is.null(schemeFastas)){

    cat(paste0('Downloading ',org,
               ' scheme ',scheme,
               ' MLST sequences at ',ddir,'/ .\n'))
    schemeFastas <- downloadPubmlst_seq(org = org,
                                        scheme=scheme,
                                        dir = ddir,
                                        n_threads = 1L)
  }else{

    if(any(!file.exists(schemeFastas))){
      warning('One or more schemeFastas do/es not exists in specified path. Downloading..',
              immediate. = T)
      cat(paste0('Downloading ',org,
                 ' scheme ',scheme,
                 ' MLST sequences at ',ddir,'/ .\n'))
      schemeFastas <- downloadPubmlst_seq(org = org,
                                          scheme=scheme,
                                          dir = ddir,
                                          n_threads = 1L)
    }else{
      all(grepl('.fas$',schemeFastas)) -> ext
      if (ext){
        for(i in 1:length(schemeFastas)){
          readLines(schemeFastas[i])-> rl
          rl[grep('>',rl)]->rl
          all(sapply(strsplit(rl,'_'),length)>=2) -> three
          if(!three){
            stop("Sequence headers doesn't have the expected format. Please download them directly
               from pubmlst.org or use the supplied R functions in this package.")
          }
        }
      }else{
        stop('Sequences should have the ".fas" extension. Please download them directly from
           pubmlst.org or use the supplied R functions in this package.')
      }
    }
  }


  #Scheme profile file: If NULL, download; else, check if ok.
  if (is.null(schemeProfile)){
    cat(paste0('Downloading ',org,
               ' scheme ',scheme,
               ' MLST profile at ',ddir,'/ .\n'))
    schemeProfile <- downloadPubmlst_profile(org = org,
                                             scheme = scheme,
                                             dir = ddir)
    prof <- read.csv(schemeProfile,sep = '\t',header = T)

  }else{
    if(!file.exists(schemeProfile)){
      warning('schemeProfile file missing, downloading..',immediate. = T)
      cat(paste0('Downloading ',org,
                 ' scheme ',scheme,
                 ' MLST profile at ',ddir,'/ .\n'))
      schemeProfile <- downloadPubmlst_profile(org = org,
                                               scheme = scheme,
                                               dir = ddir)
      prof <- read.csv(schemeProfile,sep = '\t',header = T)
    }else{
      prof <- read.csv(schemeProfile,sep = '\t',header = T)
      fils <- sub('.fas$','',sapply(schemeFastas,function(x){rev(strsplit(x,'/')[[1]])[1]}))
      if (!all(fils%in%colnames(prof))){
        stop("The scheme fasta files doesn't correspond to the supplied profile file.
      (One or more file names are not represented in the profile's column names)")
      }
    }
  }


  if (!is.null(schemeFastas)){
    #check if blast databases exists
    dds <- lapply(sub('fas$', '', schemeFastas), function(x){
      paste0(x, c('nsq', 'nin', 'nhr'))
      })
    dds<- unlist(dds)

    if (!all(file.exists(dds))){
      #Make BLAST DATABASEs
      cat('Making BLAST databases...')
      parallel::mclapply(schemeFastas,function(x){
        makeblastdb(infile = x)
      },mc.preschedule = F,mc.cores = n_threads) -> dbs
      unlist(dbs) -> dbs
      cat(' DONE!\n')
    }else{
      dbs <- sub('fas$', '', schemeFastas)
    }
  }else{
    #Make BLAST DATABASEs
    cat('Making BLAST databases...')
    parallel::mclapply(schemeFastas,function(x){
      makeblastdb(infile = x)
    },mc.preschedule = F,mc.cores = n_threads) -> dbs
    unlist(dbs) -> dbs
    cat(' DONE!\n')
  }

  #Determine mlst for each genome
  normalizePath(infiles) -> infiles
  cat('Running BLASTN...')
  resu <- parallel::mclapply(infiles,function(x){
    mlst(genome = x,
         dbs = dbs,
         write = write,
         prefix = fdir,
         dir = getwd(),
         # dnw = dnw,
         n_threads = 1L,
         outf = tempdir(),
         pid = pid,
         scov = scov)
  }, mc.preschedule = T,mc.cores = n_threads)
  cat(' DONE!\n')
  resu <- do.call(rbind,resu)
  rownames(resu) <- sub('[.]\\w+$','',basename(infiles))
  resu <- as.data.frame(resu, stringsAsFactors = FALSE)

  resu <- processResu(resu = resu,
                      write = write,
                      dbs = dbs,
                      dir = paste0(getwd(),'/'),
                      prefix = fdir)

  #Detect ST
  apply(resu,1,function(x){
    prof$ST[which(apply(prof[,colnames(resu)],1,function(y){
      all(x==y)
    }))]
  }) -> ig
  ig[!sapply(ig,length)] <- NA
  if(length(ig)>0){
    resu$ST <- unlist(ig)
  }else{
    resu$ST <- NA
  }


  #Return
  prof <- prof[, colnames(resu)]
  out <- list(result = resu,
              profile = prof)

  attr(out, 'infiles') <- basename(infiles)
  attr(out, 'org') <- org
  attr(out, 'scheme') <- scheme
  attr(out, 'write') <- write
  attr(out, 'pid') <- pid
  attr(out, 'scov') <- scov
  attr(out, 'class') <- 'mlst'

  return(out)
}


