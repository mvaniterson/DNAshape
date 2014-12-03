
##Parker figure 1G

##generate all kmers sequences for given alphabet
generateAllKmers <- function(k, alphabet=c("A", "C", "G", "T"))
  {
    require(Biostrings)
    L <- length(alphabet)
    kmers <- matrix(NA, nrow=L^k, ncol=k)
    for(i in 1:k)
      kmers[,i] <- rep(alphabet, each=L^(i-1))
    DNAStringSet(apply(kmers, 1, paste0, collapse=""))
  }

generateAllKmers <- function(k, alphabet=c("A", "C", "G", "T"))
  {
    require(Biostrings)
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

kmer7 <- generateAllKmers(7, alphabet=c("A", "C", "G", "T"))
kmer7
kmer11 <- generateAllKmers(11, alphabet=c("A", "C", "G", "T"))
kmer11

writeXStringSet(kmer7, file="kmer7seq.txt")

ifile <- file.path(getwd(), "kmer7seq.txt")
predictDNAShapes(ifile)

ofiles <- list.files(path=dirname(ifile), pattern="kmer7seq")[-1]
ofiles
i <- 2

##too slow use linux
##shape <- readShape(file.path(dirname(ifile), ofiles[i]))

file <- file.path(dirname(ifile), ofiles[i])
shape <- read.table(pipe(paste("sed '/>/d'", file)), sep=",")
dim(shape)

str(shape)
dim(shape)
head(shape)

reducedShape <- shape[, c(-c(1:2), -c((ncol(shape)-1):ncol(shape)))]
dim(reducedShape)

##figure 1a,b,c
library(DNAshape)
library(Biostrings)

seqsA <- DNAStringSet(c("ATACGCG", "ATAGGCG", "ATATGCG"))
names(seqsA) <- paste0("sequence", 1:3)
plotParker(seqsA)

seqsB <- DNAStringSet(c("ACGTACC", "AGGGAGC", "ATGCATC"))
names(seqsB) <- paste0("sequence", 1:3)
plotParker(seqsB)

seqsC <- DNAStringSet(c("AATGTTT", "AATTTTT", "AATCTTT"))
names(seqsC) <- paste0("sequence", 1:3)
plotParker(seqsC)

plotParker<- function(seqs)
  {
    faFile <- "Parkerseqs.fa"
    writeXStringSet(DNAStringSet(seqs), file=faFile)
    predictDNAShapes(file.path(getwd(), faFile))

    ##plot single shape
    shapeFiles <- list.files(path=dirname(faFile), pattern="seq")[-1]
    shapeFiles
    op <- par(mfcol=c(2,2))
    for(i in 1:4)
      {
        shape <- readShape(file.path(dirname(faFile), shapeFiles[i]))
        shapeType <- gsub(paste0(faFile, "."), "", basename(shapeFiles[i]))
        plotShape(shape, seqs, shapeType=shapeType)
      }
    par(op)
  }
