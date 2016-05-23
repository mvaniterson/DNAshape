##' generate all kmers of given length
##'
##' generate all kmers of given length
##' with middle nucleotide...
##' TODO find a fast way to exclude reverse complement sequences
##' @title generate kmers
##' @param k length
##' @param alphabet default DNA bases 
##' @return DNAStringSet with all kmers 
##' @author mvaniterson
##' @export
##' @importFrom Biostrings DNAStringSet
kmers <- function(k, alphabet=c("A", "C", "G", "T")) {
    L <- length(alphabet)
    kmers <- matrix(NA, nrow=L^k, ncol=k)
    for(i in 1:k)
      kmers[,i] <- rep(alphabet, each=L^(i-1))

    ##swap first column to the middle
    col1 <- kmers[,1]
    col2 <- kmers[,floor(k/2)+1]

    kmers[,1] <- col2
    kmers[,floor(k/2)+1] <- col1

    DNAStringSet(apply(kmers, 1, paste0, collapse=""))
  }
