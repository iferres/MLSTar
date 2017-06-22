# quickMLST

This R package allows you to easily determine the Multi Locus Sequence Type (MLST) of your genomes. It also works as an interface between [PubMLST](https://pubmlst.org/) through their [RESTful API](https://pubmlst.org/rest/), so you don't have to bother downloading and collecting files: the application does it automatically.

NOTE: A good internet connection is highly recommended when using this package. A delay or a malfunction may be due to problems with the pubmlst.org server. 

## External dependencies



## Quick standard workflow

The first step in your analysis should be to check the in pubmlst.org database if your organism of interest is available. So, first load the package and then run `listPubmlst_orgs()` function, printing only the first 50 elements:

``` r
library(quickMLST)
listPubmlst_orgs()[1:50]
```

    ##  [1] "achromobacter"           "abaumannii"             
    ##  [3] "aeromonas"               "aphagocytophilum"       
    ##  [5] "arcobacter"              "bcereus"                
    ##  [7] "blicheniformis"          "bsubtilis"              
    ##  [9] "bhenselae"               "bordetella"             
    ## [11] "borrelia"                "brachyspira"            
    ## [13] "brucella"                "bcc"                    
    ## [15] "bpseudomallei"           "campylobacter_nonjejuni"
    ## [17] "campylobacter"           "calbicans"              
    ## [19] "cglabrata"               "ctropicalis"            
    ## [21] "cmaltaromaticum"         "chlamydiales"           
    ## [23] "cfreundii"               "csinensis"              
    ## [25] "cbotulinum"              "cdifficile"             
    ## [27] "csepticum"               "cdiphtheriae"           
    ## [29] "cronobacter"             "dnodosus"               
    ## [31] "ecloacae"                "efaecalis"              
    ## [33] "efaecium"                "fpsychrophilum"         
    ## [35] "hinfluenzae"             "hparasuis"              
    ## [37] "hcinaedi"                "helicobacter"           
    ## [39] "hsuis"                   "koxytoca"               
    ## [41] "kseptempunctata"         "lsalivarius"            
    ## [43] "leptospira"              "mcanis"                 
    ## [45] "mcaseolyticus"           "mplutonius"             
    ## [47] "mabscessus"              "mycobacteria"           
    ## [49] "magalactiae"             "mbovis"

Lets say we are interested in Leptospira genus, which is in the place 43 in the list above. So:

``` r
listPubmlst_orgs() -> lst
lst[43]
```

    ## [1] "leptospira"

Now, lets check for available MLST schemes for this organism:

``` r
listPubmlst_schemes(org = lst[43])
```

    ## $scheme_1
    ## [1] "glmU_1" "pntA_1" "sucA_1" "tpiA_1" "pfkB_1" "mreA_1" "caiB_1"
    ## attr(,"Desc")
    ## [1] "MLST (scheme 1)"
    ## 
    ## $scheme_2
    ## [1] "adk_2"    "glmU_2"   "icdA_2"   "lipL32_2" "lipL41_2" "mreA_2"  
    ## [7] "pntA_2"  
    ## attr(,"Desc")
    ## [1] "MLST (scheme 2)"
    ## 
    ## $scheme_3
    ## [1] "adk_3"    "icdA_3"   "lipL32_3" "lipL41_3" "rrs2_3"   "secY_3"  
    ## attr(,"Desc")
    ## [1] "MLST (scheme 3)"

As you can see, `listPubmlst_schemes` return a list with the loci names corresponding to each scheme. As an attribute of each list element there is information about each mlst scheme.

Now you can choose between two ways: the easy way and the hard way.

The hard way implies calling `downloadPubmlst_seq(org = lst[43], scheme = 1)` and then `downloadPubmlst_profile(org = lst[43], scheme = 1)` functions included in this package to download the scheme fasta files and the profile tab file for the organism and the scheme of interest, and then passing the files to the subsequent `doMLST()` function to `schemeFastas` and `schemeProfile` arguments.

The easy way is to left those arguments `NULL` (default), and let the `doMLST()` function do it for you.

Let see an example with toy data attached on this package:

``` r
#First we list the atteched tar.gz file
system.file('extdata', 'toyExample.tar.gz', package = 'quickMLST') -> tgz
untar(tarfile = tgz, exdir = getwd(), list = T) -> genomes
#Decompress them
untar(tarfile = tgz,exdir = getwd())
genomes
```

    ## [1] "L._borgpetersenii.fna" "L._interrogans.fna"    "L._kirschneri.fna"

In this example we have 3 pathogenic leptospira genomes, in fasta format.

Lets determine the MLST for the scheme 3.

``` r
doMLST(infiles = genomes, # The fasta files
       org = lst[43], # The organism, in this case is "leptospira"
       scheme = 3, # Scheme id number
       write.new = FALSE, # Don't write fasta files for new alleles found
       dir = getwd(), # Put MLST allele files in this dir
       n_threads = 3) -> res # Use 3 threads
```

    ## Downloading leptospira scheme 3 MLST sequences at /home/iferres/Documents/mlst_Lepto// .
    ## Downloading leptospira scheme 3 MLST profile at /home/iferres/Documents/mlst_Lepto// .
    ## Making BLAST databases... DONE!
    ## Running BLASTN... DONE!

``` r
#Output:
res
```

    ##                       adk_3 icdA_3 lipL32_3 lipL41_3 rrs2_3 secY_3 ST
    ## L._borgpetersenii.fna    57     54        u       39     20     47 NA
    ## L._interrogans.fna        2      2        2        2      1      2 47
    ## L._kirschneri.fna        13     24       11       16     12     22 94

As you can see, a `data.frame` is returned. Each row is a genome, and each column is a scheme locus. The number refers to the allele number id.

A `"u"` means that a new allele was found, e.g. `res$lip32_3[1]`: this allele is not yet reported in the pubmlst database. If option `write.new` is set to `TRUE`, then a fasta file is written in `dir` with this new allele.

A `<NA>` means that no allele was found, i.e. no blastn local alignment pass the inclusion threshold (by default, this threshold are a percentage identity grater or equal to 90, and a subject coverage greater or equal to 0.9). In this example this was no the case for any of the screened genomes.

The last column refers to the Sequence Type (ST). If possible, the function identifies the ST of each genome, otherwise a `NA` is returned (e.g. `res$ST[1]`).

An easy way of obtaining the composition of the 3 mlst schemes available for this organism would be:

``` r
lapply(1:3,function(x){
  
  doMLST(infiles = genomes, # The fasta files
         org = lst[43], # The organism, in this case is "leptospira"
         scheme = x, # Scheme id number. Will iterate between 1 and 3.
         write.new = FALSE, # Don't write fasta files for new alleles found
         dir = getwd(), # Put MLST allele files in this dir
         n_threads = 3)
  
}) -> allres
```

    ## Downloading leptospira scheme 1 MLST sequences at /home/iferres/Documentos/mlst_Lepto// .
    ## Downloading leptospira scheme 1 MLST profile at /home/iferres/Documentos/mlst_Lepto// .
    ## Making BLAST databases... DONE!
    ## Running BLASTN... DONE!
    ## Downloading leptospira scheme 2 MLST sequences at /home/iferres/Documentos/mlst_Lepto// .
    ## Downloading leptospira scheme 2 MLST profile at /home/iferres/Documentos/mlst_Lepto// .
    ## Making BLAST databases... DONE!
    ## Running BLASTN... DONE!
    ## Downloading leptospira scheme 3 MLST sequences at /home/iferres/Documentos/mlst_Lepto// .
    ## Downloading leptospira scheme 3 MLST profile at /home/iferres/Documentos/mlst_Lepto// .
    ## Making BLAST databases... DONE!
    ## Running BLASTN... DONE!

``` r
allres
```

    ## [[1]]
    ##                       glmU_1 pntA_1 sucA_1 tpiA_1 pfkB_1 mreA_1 caiB_1  ST
    ## L._borgpetersenii.fna     26     30     28     35     39     29     29 152
    ## L._interrogans.fna         1      1      1      1      1      1      1   1
    ## L._kirschneri.fna         19     20     13     22     31     18     23 110
    ## 
    ## [[2]]
    ##                       adk_2 glmU_2 icdA_2 lipL32_2 lipL41_2 mreA_2 pntA_2
    ## L._borgpetersenii.fna    29      u      u        u        u     31      u
    ## L._interrogans.fna        3      1      2        2        5      5      2
    ## L._kirschneri.fna        28     25     21        8        7      7     11
    ##                        ST
    ## L._borgpetersenii.fna  NA
    ## L._interrogans.fna      7
    ## L._kirschneri.fna     100
    ## 
    ## [[3]]
    ##                       adk_3 icdA_3 lipL32_3 lipL41_3 rrs2_3 secY_3 ST
    ## L._borgpetersenii.fna    57     54        u       39     20     47 NA
    ## L._interrogans.fna        2      2        2        2      1      2 47
    ## L._kirschneri.fna        13     24       11       16     12     22 94

That's it. Now we have the MLST of our genomes for the 3 available schemes.

You should check the files downloaded from [PubMLST](https://pubmlst.org/) on your working directory .

## Installation

The easiest way to install the package is from within R, using `devtools`:

```r
library(devtools)
install_github('iferres/quickMLST')
```

Alternatively, you can clone the repository and install it manually:

Change directory to where you want to download the package (e.g. ~/Downloads/), and clone the repository to your machine:

```
cd ~/Downloads/
git clone https://github.com/iferres/quickMLST
```
Then, from R, install the package:

```r
install.packages(pkgs = '~/Downloads/quickMLST', type = 'source', repos = NULL)
```

### External dependencies

This package uses [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) as search engine, so it should be installed prior to running the functions.

