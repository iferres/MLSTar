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
#' @param prefix \code{character} A prefix to the fasta files of found
#' sequences (see \code{write}). (Default: \code{"allele_"}).
#' @param dir An existing directory where to put the loci fasta files in case
#' they are not provided by the user. Also sequences found will be placed here
#' if \code{write} is set to ethier \code{"new"} or \code{"all"}.
#' @param n_threads \code{integer}. The number of cores to use. Each job consist
#' on the process for one genome. Blastn searches will use 1 core per job.
#' @param outf Where blastn output will be written. By default, in a temp
#' directory.
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
#' @importFrom parallel mclapply
#' @export
doMLST <- function(infiles,
                   org='leptospira',
                   scheme=1L,
                   schemeFastas=NULL,
                   schemeProfile=NULL,
                   write = 'new',
                   prefix = 'allele_',
                   dir='.',
                   n_threads=1L,
                   outf=tempdir(),
                   pid=90L,
                   scov=0.9){

  if (.Platform$OS.type!='unix'){
    stop('Sorry, this package works only on unix-like platforms.')
  }

  if(Sys.which('blastn')==''){
    stop('blastn binary is not in $PATH. Please install it before running this function.')
  }

  paste0(normalizePath(dir),'/') -> dir

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
               ' MLST sequences at ',dir,'/ .\n'))
    schemeFastas <- downloadPubmlst_seq(org = org,
                                        scheme=scheme,
                                        dir = dir,
                                        n_threads = 1L)
  }else{

    if(any(!file.exists(schemeFastas))){
      warning('One or more schemeFastas do/es not exists in specified path. Downloading..',
              immediate. = T)
      cat(paste0('Downloading ',org,
                 ' scheme ',scheme,
                 ' MLST sequences at ',dir,'/ .\n'))
      schemeFastas <- downloadPubmlst_seq(org = org,
                                          scheme=scheme,
                                          dir = dir,
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
               ' MLST profile at ',dir,'/ .\n'))
    schemeProfile <- downloadPubmlst_profile(org = org,
                                             scheme = scheme,
                                             dir = dir)
    prof <- read.csv(schemeProfile,sep = '\t',header = T)

  }else{
    if(!file.exists(schemeProfile)){
      warning('schemeProfile file missing, downloading..',immediate. = T)
      cat(paste0('Downloading ',org,
                 ' scheme ',scheme,
                 ' MLST profile at ',dir,'/ .\n'))
      schemeProfile <- downloadPubmlst_profile(org = org,
                                               scheme = scheme,
                                               dir = dir)
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


  #Make BLAST DATABASEs
  cat('Making BLAST databases...')
  parallel::mclapply(schemeFastas,function(x){
    makeblastdb(infile = x)
  },mc.preschedule = F,mc.cores = n_threads) -> dbs
  unlist(dbs) -> dbs
  cat(' DONE!\n')


  #Determine mlst for each genome
  normalizePath(infiles) -> infiles
  cat('Running BLASTN...')
  parallel::mclapply(infiles,function(x){
    mlst(genome = x,
         dbs = dbs,
         write = write,
         prefix = prefix,
         dir = dir,
         n_threads = 1L,
         outf = outf,
         pid = pid,
         scov = scov)
  }, mc.preschedule = T,mc.cores = n_threads) -> ress
  cat(' DONE!\n')
  do.call(rbind,ress) -> resu
  rownames(resu) <- sub('.fas$','',sapply(infiles,function(x){rev(strsplit(x,'/')[[1]])[1]}))
  as.data.frame(resu) -> resu

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


  #Remove blast databases
  c('.nsq','.nin','.nhr')-> dbext
  unlist(sapply(sub('.fas$','',schemeFastas),function(x){
    paste0(x,dbext)
    },simplify = F)) -> el
  file.remove(el)

  #Return
  return(resu)
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
#' sequences (see \code{write}). (Default: \code{"allele_"}).
#' @param dir Sequences of new allele found will be placed here.
#' @param n_threads \code{integer}. The number of threads to use by BLASTN.
#' @param outf Where the blastn output will be written.
#' @param pid Percentage identity.
#' @param scov Query coverage.
#' @return A vector with the mlst.
#' @importFrom parallel mclapply
#' @importFrom utils read.csv
#' @author Ignacio Ferres
mlst <- function(genome,
                 dbs,
                 write='new',
                 prefix = 'allele_',
                 dir='.',
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
                            pid = pid,
                            scov = scov)
    res[i] <- as.character(a)
    nams[i] <- names(a)
  }

  names(res) <- nams
  res
}




