.seq2shape  <- function(sequence, kmerLength=5){
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
    shape <- shape.profiles[match(kmers, rownames(shape.profiles)), properties, drop=FALSE]

    ##this is not really nice
    ##average pair properties
    z <- c(NA, NA)
    for(i in 3:4) {
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



##' predict shape for DNA sequence
##'
##' R implementation that query's the extended
##' (reverse complement) pentamer shape profiles
##' TODO add fasta-file input
##' @title seq2shape
##' @param sequence character DNA sequence
##' (reverse complement) pentamer shape profiles
##' @param kmerLength default = 5
##' @return shape profile
##' @author mvaniterson
##' @importFrom graphics abline axTicks mtext par plot points
##' @importFrom utils data read.table
##' @export
seq2shape  <- function(sequence, kmerLength=5) {
    if(class(sequence) == "character" & length(sequence) == 1)
        return(.seq2shape(sequence))
    else if(class(sequence) == "character" & length(sequence) > 1)
        return(lapply(sequence, .seq2shape))
    else if(class(sequence) == "DNAString")
        return(.seq2shape(as.character(sequence)))
    else if(class(sequence) == "DNAStringSet")
        return(lapply(as.character(sequence), .seq2shape))
    else
        stop(paste("Unknown sequence input:", class(sequence)))
}

##' calculate shape difference
##' 
##' Calculate shape difference between two sequences
##' @title shapediff
##' @param seq1 DNA sequence 1
##' @param seq2 DNA sequence 2
##' @param type shape difference measure rmsd or eucl.
##' @return shape difference
##' @author mvaniterson
##' @export
shapediff <- function(seq1, seq2, type=c("rmsd", "eucl")) {
    type <- match.arg(type)    
    if(nchar(as.character(seq1)) != nchar(as.character(seq2)))
        stop("Sequence lengths do not match!")

    shape1 <- seq2shape(seq1)
    shape2 <- seq2shape(seq2)

    switch(type,
           eucl = colSums((shape1 - shape2)^2, na.rm=TRUE),
           rmsd = sqrt(colSums((shape1 - shape2)^2, na.rm=TRUE)/nrow(shape1)))
}
