    ############################
#### INTERFACE WITH PUBMLST.ORG ####
    ############################

#' @name listPubmlst_orgs
#' @title List Available Organisms in Pubmlst.org Database
#' @description List all organisms available in pubmlst database.
#' @param try.times \code{integer}. Try n times to connect with pubmlst.org before returning an error message.
#' @return A list of genus.
#' @importFrom httr GET status_code
#' @importFrom jsonlite fromJSON
#' @author Ignacio Ferres
#' @references Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming
#' interface for the PubMLST molecular typing and genome databases. Database:
#' The Journal of Biological Databases and Curation. 2017;2017:bax060.
#' doi:10.1093/database/bax060.
#' @export
listPubmlst_orgs <- function(try.times=3){
  cnt <- 0
  status.code <- 0
  while (cnt<try.times && status.code!=200){
    httr::GET('http://rest.pubmlst.org/') -> g
    httr::status_code(g) -> status.code
    cnt <- cnt + 1
  }
  if (httr::status_code(g)==200){

    rawToChar(g$content) -> a
    jsonlite::fromJSON(a)$databases -> b
    lapply(b,function(x){
      sub('_seqdef','',sub('pubmlst_','',x$name[which(grepl('seqdef',x$name))]))
    }) -> b
    unlist(b)

  }else{
    stop('Something went wrong with the connection.')
  }
}

#' @name listPubmlst_schemes
#' @title List Available MLST Schemes for a Given Genus
#' @description Takes a single genus as returned by \link{listPubmlst_orgs},
#' and lists the available mlst schemes for that genus in pubmlst.org .
#' @param org \code{character}. A genus (lowercase).
#' @return A \code{list} with genes for each available mlst scheme.
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' @author Ignacio Ferres
#' @references Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming
#' interface for the PubMLST molecular typing and genome databases. Database:
#' The Journal of Biological Databases and Curation. 2017;2017:bax060.
#' doi:10.1093/database/bax060.
#' @export
listPubmlst_schemes <- function(org='leptospira'){
  paste0('http://rest.pubmlst.org/db/pubmlst_',org,'_seqdef/schemes') -> d
  jsonlite::fromJSON(rawToChar(httr::GET(d)$content))$schemes -> e
  if(is.null(e$scheme)){
    stop('Invalid organism.')
  }
  sapply(e$scheme,function(x){rev(strsplit(x,'/')[[1]])[1]}) -> nscheme
  lapply(e$scheme,function(x){jsonlite::fromJSON(rawToChar(httr::GET(x)$content))$loci}) -> lc
  lapply(lc,function(x){do.call(rbind,strsplit(x,'/'))[,7]}) -> gns
  names(gns) <- paste0('scheme_',nscheme)
  for (i in 1:length(gns)){
    attr(gns[[i]],'Desc') <- e$description[i]
  }
  gns
}

#' @name downloadPubmlst_seq
#' @title Download Genes for a Given MLST Scheme
#' @description This function takes a genus, as returned by \link{listPubmlst_orgs},
#' and a mlst scheme as either an \code{integer} or as a vector of locus names
#' as returned by \link{listPubmlst_schemes}, and downloads the fasta files for
#' the given genus and mlst scheme.
#' @param org \code{character}. A genus (lowercase).
#' @param scheme Either an \code{integer} specifiying the scheme number to
#' download, or a vector of the locus names of that organism.
#' @param dir Where to put the fasta files.
#' @param n_threads \code{integer}. The number of cores to use.
#' @return A vector of paths to the downloaded fasta files.
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' @importFrom parallel mclapply
#' @author Ignacio Ferres
#' @references Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming
#' interface for the PubMLST molecular typing and genome databases. Database:
#' The Journal of Biological Databases and Curation. 2017;2017:bax060.
#' doi:10.1093/database/bax060.
#' @importFrom utils write.table
#' @export
downloadPubmlst_seq <- function(org='leptospira',
                                scheme='',
                                dir='.',
                                n_threads=1L){

  paste0(normalizePath(dir),'/') -> dir

  if(any(class(scheme)==c('numeric','integer'))){
    paste0('http://rest.pubmlst.org/db/pubmlst_',org,'_seqdef/schemes/',scheme) -> d
    jsonlite::fromJSON(rawToChar(httr::GET(d)$content)) -> jsn
    jsn$loci -> lc
    if (is.null(lc)){
      stop('Could not download sequences - Invalid input.')
    }
    do.call(rbind,strsplit(lc,'/'))[,7] -> scheme
  }

  paste0('http://rest.pubmlst.org/db/pubmlst_',org,'_seqdef/loci/',scheme,'/alleles_fasta') -> d

  parallel::mclapply(d,function(x){
    strsplit(rawToChar(httr::GET(x)$content),'\n')[[1]]
  },mc.preschedule = F,mc.cores = n_threads) -> lp
  paste0(dir,scheme,'.fas')-> nam
  for (i in 1:length(lp)){
    seqinr::write.fasta(sequences = as.list(lp[[i]][seq(2,length(lp[[i]]),2)]),
                        names = as.list(sub('>','',lp[[i]][seq(1,length(lp[[i]]),2)])),
                        file.out = nam[i],
                        as.string = T)
  }

  nam
}

#' @name downloadPubmlst_profile
#' @title Download the Profile for a Given MLST Scheme
#' @description This function takes a genus, as returned by \link{listPubmlst_orgs},
#' and a mlst scheme as an \code{integer} and downloads the profile file for
#' the given genus and mlst scheme, in tsv format.
#' @param org \code{character}. A genus (lowercase).
#' @param scheme An \code{integer} specifiying the scheme number.
#' @param dir Where to put the tsv profile file.
#' @return A vector of paths to the downloaded tsv file.
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' @author Ignacio Ferres
#' @references Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming
#' interface for the PubMLST molecular typing and genome databases. Database:
#' The Journal of Biological Databases and Curation. 2017;2017:bax060.
#' doi:10.1093/database/bax060.
#' @export
downloadPubmlst_profile <- function(org='leptospira',
                                    scheme=1L,
                                    dir='.'){

  paste0(normalizePath(dir),'/') -> dir
  paste0('http://rest.pubmlst.org/db/pubmlst_',org,'_seqdef/schemes/',scheme) -> d
  jsonlite::fromJSON(rawToChar(httr::GET(d)$content))$profiles_csv -> d2
  if (is.null(d2)){
    stop('Could not download profile - Invalid input.')
  }
  rawToChar(httr::GET(d2)$content) -> f
  paste0(dir,'profile_scheme',scheme,'.tab') -> out
  write.table(f,
              file = out,
              quote = F,
              sep = '\t',
              row.names = F,
              col.names = F)
  out

}
