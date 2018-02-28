#Internals

    ########
#### BLASTN ####
    ########

#' @name makeblastdb
#' @title Make a Blast Database
#' @description Takes a nucleotide fasta file and produce a blast database.
#' @param infile A nucleotide fasta file.
#' @param hash_index \code{logical}. See Blast+ manual.
#' @param parse_seqids \code{logical}. See Blast+ manual.
#' @return The path to the created database.
#' @references  Altschul, Gish, Miller, Myers & Lipman. (1990) "Basic local
#'  alignment search tool." \emph{J. Mol. Biol}. \strong{215}:403-410.
#' @author Ignacio Ferres
makeblastdb <- function(infile,
                        hash_index=F,
                        parse_seqids=F){
  args <- c('-dbtype nucl')
  if(hash_index){
    args <- paste(args,'-hash_index')
  }
  if(parse_seqids){
    args <- paste(args,'-parse_seqids')
  }

  nam <- sub('[.]\\w+$','',infile)
  cmd<-paste('makeblastdb -in',infile,
             args,
             '-out',nam,
             '-title',nam)
  system(cmd,ignore.stdout = TRUE)
  nam

}


#' @name blastn
#' @title Seach with Blastn.
#' @description Seach nucleotide sequences in a nucleotide database.
#' @param genome The query. A nucleotide fasta file.
#' @param db The subject. A nucleotide blast database as created by
#' \link{makeblastdb}.
#' @param outdir Where to put the blast output file.
#' @param n_threads \code{integer}. The nomber of cores to use.
#' @param eval Evalue reporting threshold.
#' @return The path to the blastn output. The blastn output is written in the
#' following format as specified by \code{-outfmt} option in Blast+ program:
#'
#' \code{-outfmt '6 qseqid sseqid pident gaps length slen evalue qseq'}.
#' @references  Altschul, Gish, Miller, Myers & Lipman. (1990) "Basic local
#'  alignment search tool." \emph{J. Mol. Biol}. \strong{215}:403-410.
#' @author Ignacio Ferres
blastn <- function(genome='',
                   db='',
                   outdir='.',
                   n_threads=1L,
                   eval=1e-6){

  paste0(normalizePath(outdir),'/') -> outdir
  dbn <- rev(strsplit(db,'/')[[1]])[1]
  gnam <- rev(strsplit(genome,'/')[[1]])[1]
  outfile <- paste0(outdir,sub('.f\\w+$','',gnam),'_vs_',sub('[.]\\w+$','',dbn))
  cmd<-paste0("blastn -word_size 50 -ungapped -dust no -query ",genome,
              " -db ",db,
              " -evalue ",eval,
              " -outfmt '6 qseqid sseqid pident gaps length qstart qend evalue qseq'",
              " -out ",outfile,
              " -num_threads ",n_threads)
  system(cmd)
  outfile

}


#' @name readBlastResult
#' @title Read Blast Result
#' @description Read blast result
#' @param blout The path to the blastn output, written in format as specified
#' by the Blast+ parameter \code{-outfmt '6 qseqid sseqid pident gaps length
#' slen evalue qseq'}.
#' @return A \code{data.frame} or a \code{try-error} message if blout is empthy.
#' @author Ignacio Ferres
#' @importFrom utils read.table
readBlastResult <- function(blout=''){

  cols<-c('qid','sid','pid','gaps','lgth','qstart','qend','evalue','qseq')
  try(read.csv(blout,
               header = F,
               col.names = cols,
               sep = '\t',
               stringsAsFactors = FALSE),
      silent = T) -> res
  rev(strsplit(blout,'/')[[1]])[1] -> fi
  strsplit(fi,'_vs_')[[1]][1] -> infile
  attr(res,'infile') <- infile
  strsplit(fi,'_vs_')[[1]][2] -> indb
  attr(res,'indb') <- indb
  res

}


#' @name processBlastResult
#' @title Process Blastn Result
#' @description Process blastn output. Reports whether the genome
#' contains a reported allele in pubmlst, or it has a not reported (new)
#' allele, or if no allele where found. In case a new allele were found, a
#' fasta file is written.
#' @param blastRes A \code{data.frame} as returned by \link{readBlastResult}.
#' @param pid Percentage identity reporting threshold.
#' @param scov Query coverage reporting threshold.
#' @param write \code{character}. One of \code{"new"} (Default), \code{"all"} or
#' \code{"none"}. The fist one writes only new alleles found (not reported in
#' pubmlst.org), the second writes all alleles found, and "none" do not write
#' any file.
#' @param prefix \code{character} A prefix to the fasta files of found
#' sequences (see \code{write}). (Default: \code{"allele"}).
#' @param dir The directory where to put the fasta file of allele sequences
#' found.
#' @param dnw Temporary directory where to put new alleles to check later if
#' are the same.
#' @return The allele number id if an exact match with a reported mlst gene is
#' found, a name arbitrarily given if a new allele is found, or \code{NA} if no
#' alleles are found.
#' @importFrom seqinr write.fasta
#' @author Ignacio Ferres
processBlastResult <- function(blastRes,
                               pid=90,
                               scov=0.9,
                               write='new',
                               prefix = 'alleles',
                               dir='.',
                               dnw=paste0(dir, 'tmp/')
                               ){

  write <- match.arg(write, c('none', 'new', 'all'))
  dir <- paste0(normalizePath(dir),'/')
  gene <- attr(blastRes,'indb')
  out.newAllele <- paste0(dir, prefix, '.', gene, '.fasta')
  out.tmp <- paste0(dnw, prefix, '.', gene, '.fasta')
  gid <- attr(blastRes,'infile')

  if (class(blastRes)!='try-error'){
    blastRes$scov <- (blastRes$lgth - blastRes$gaps) /
      (blastRes$qend - blastRes$qstart + 1)

    if(any(blastRes$pid==100 &
           blastRes$scov==1)){

      wh <- which(blastRes$pid==100 &
                    blastRes$scov==1)
      hit <- blastRes$sid[wh]
      spl <- rev(strsplit(hit,'_')[[1]])[1]
      allele <- as.character(spl)
      names(allele) <- gene

      if(write=='all'){
        sq <- blastRes$qseq[wh]
        qid <- blastRes$qid[wh]
        nsq <- paste0(hit, ';', gid, ';', qid)
        write.fasta(sequences = sq,
                    names = nsq,
                    file.out = out.newAllele,
                    open = 'a',
                    as.string = TRUE)
      }

      return(allele)

    }else if(any(blastRes$pid>=pid &
                 blastRes$scov>=scov)){

      blastRes$Val <- blastRes$pid/blastRes$scov
      wh <- which.min(abs(blastRes$Val - 100))

      allele <- 'u'
      names(allele) <- gene

      sq <- blastRes$qseq[wh]
      qid <- blastRes$qid[wh]
      hit <- paste0(gene, '_u')
      nsq <- paste0(hit, ';', gid, ';', qid)

      if (write%in%c('new', 'all')){

        write.fasta(sequences = sq,
                    names = nsq,
                    file.out = out.newAllele,
                    open = 'a',
                    as.string = TRUE)

      }

      # Write anyway to dir/tmp/ to check later if are the same
      write.fasta(sequences = sq,
                  names = nsq,
                  file.out = out.tmp,
                  open = 'a',
                  as.string = TRUE)


      return(allele)

    }else{

      allele <- NA
      names(allele) <- gene
      return(allele)

    }


  }else{

    allele <- NA
    names(allele) <- gene
    return(allele)

  }


}
