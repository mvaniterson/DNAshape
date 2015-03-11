##' predict shape for DNA sequence
##'
##' R implementation that query's the extended
##' (reverse complement) pentamer shape profiles
##' @title seq2shape
##' @param sequence character DNA sequence
##' @param kmerShapes if missing load's the extended
##' (reverse complement) pentamer shape profiles
##' @param kmerLength default = 5
##' @return shape profile
##' @author mvaniterson
seq2shape  <- function(sequence, kmerShapes=NULL, kmerLength=5)
  {
    
    if(is.null(kmerShapes))
      data(shape.profiles, package="DNAshape")
    
    ##extract kmers from sequence
    kmers <- mapply(FUN=function(x, y) substr(sequence, x, y),
                    1:(nchar(sequence)-kmerLength+1),
                    kmerLength:nchar(sequence))

    properties <- c("minor_ave",
                    "propel_ave",
                    "twist1_ave",
                    "twist2_ave",
                    "roll1_ave",
                    "roll2_ave")

    ##extract shape properties for each kmer
    shape <- kmerShapes[match(kmers, rownames(kmerShapes)), properties, drop=TRUE]

    ##this is not really nice
    ##average pair properties
    z <- c(NA, NA)
    for(i in 3:4)
      {
        x <- as.numeric(shape[,i])
        y <- as.numeric(shape[,i+1])
        shape[1, i] <- x[1]
        if(length(x) > 1)
          shape[-1,i] <- 0.5*(x[-1] + y[-length(y)])
        z <- c(z, y[length(y)])
        shape <- shape[,-i-1,drop=FALSE]
      }
    
    ##reshape output adding step bases and pairs
    shape <- rbind(shape, N=z)
    colnames(shape) <- gsub("1", "", colnames(shape))
    shape  
  }

