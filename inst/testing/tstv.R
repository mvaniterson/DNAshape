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

kmer <- generateAllKmers(7, alphabet=c("A", "C", "G", "T"))

##transitions
n <- length(kmer)
TS1 <- kmer[seq(1,n,2)] ##extract transitions ## A - G
names(TS1) <- paste0(1:length(TS1), "Y")
TS2 <- kmer[seq(2,n,2)] ##extract transitions ## C - T
names(TS2) <- paste0(1:length(TS1), "U")

##transversions
TV1 <- kmer[sort(c(seq(1, n, 4), seq(2, n, 4)))] ## A - C
names(TV1) <- paste0(1:length(TV1), c("Y","R"))
TV2 <- kmer[sort(c(seq(1, n, 4), seq(4, n, 4)))] ## A - T
names(TV2) <- paste0(1:length(TV2), c("Y","R"))
TV3 <- kmer[sort(c(seq(2, n, 4), seq(3, n, 4)))] ## G - C
names(TV3) <- paste0(1:length(TV1), c("Y","R"))
TV4 <- kmer[sort(c(seq(3, n, 4), seq(4, n, 4)))] ## G - T
names(TV4) <- paste0(1:length(TV2), c("Y","R"))

library(DNAshape)

shapeDistance <- function(seqs, tmp.dir)
  {
    writeXStringSet(seqs, file.path(tmp.dir, "tmp.fa"))
    faFile <- file.path(tmp.dir, "tmp.fa")
    predictDNAShapes(faFile)
    shapeFiles <- list.files(path=dirname(faFile), pattern="tmp")[-1]
    d <- matrix(NA, nrow=length(seqs)/2, ncol=4)
    for(i in 1:4)
      {
        shape <- readShape(file.path(dirname(faFile), shapeFiles[i]))
        shape <- do.call(rbind.data.frame, shape)
        colnames(shape) <- 1:ncol(shape)
        shape <- shape[, !is.na(shape[1,])]
        d[, i] <- sqrt(rowSums((shape[seq(1, nrow(shape), 2),] - shape[seq(2, nrow(shape), 2),])^2))  ##euclidian distiance
      }
    colnames(d) <- gsub("^.*\\.", "", shapeFiles)
    rownames(d) <- paste(seqs[seq(1, nrow(shape), 2)], seqs[seq(2, nrow(shape), 2)], sep="-")
    d
  }

dTS1 <- shapeDistance(TS1, tmp.dir ="/home/mviterson")
dTS2 <- shapeDistance(TS2, tmp.dir ="/home/mviterson")
dTV1 <- shapeDistance(TV1, tmp.dir ="/home/mviterson")
dTV2 <- shapeDistance(TV2, tmp.dir ="/home/mviterson")
dTV3 <- shapeDistance(TV3, tmp.dir ="/home/mviterson")
dTV4 <- shapeDistance(TV4, tmp.dir ="/home/mviterson")

##faster reading if sequences are of equal length
readShape <- function(file)
  {
    lines <- readLines(file)
    shapes <- lines[seq(2, length(lines), 2)]
    names(shapes) <- lines[seq(1, length(lines), 2)]
    row.names <- names(shapes)
    shapes <- t(simplify2array(lapply(shapes, function(x) unlist(strsplit(x, ",")))))
    shapes <- apply(shapes, 2, as.numeric)
    rownames(shapes) <- row.names
    print(paste("Read", nrow(shapes), "shapes."))
    shapes
  }

data <- data.frame(Mutation = factor(rep(c("A<->G", "C<->T", "A<->C", "A<->T", "G<->C", "G<->T"), each=nrow(dTS1)), levels=c("A<->G", "C<->T", "A<->C", "A<->T", "G<->C", "G<->T")),
                   MutationType = rep(c("Transition", "Transition", "Transversion", "Transversion", "Transversion", "Transversion"), each=nrow(dTS1)),
                   ShapeProfile = rep(colnames(dTS1), each=6*nrow(dTS1)),
                   Distance = c(as.vector(dTS1), as.vector(dTS2),
                     as.vector(dTV1), as.vector(dTV2),
                     as.vector(dTV3), as.vector(dTV4)))

library(ggplot2)
gp <- ggplot(data, aes(Mutation, Distance, fill=MutationType))
gp <- gp + geom_boxplot()
gp <- gp + facet_grid(ShapeProfile~., scale="free")
gp

ggsave(filename="tvtskmer9.pdf")

i <- 2
plot(dTS1[,i], dTS2[,i], main=colnames(dTV1)[i])
grid()
abline(0, 1)

writeXStringSet(TS1, file.path("/home/mviterson/", "tmp.fa"))
faFile <- file.path("/home/mviterson/", "tmp.fa")

plotShapes <- function (faFile, subset=NULL) 
{
    sequences <- readDNAStringSet(faFile, format = "fasta")
    shapeTypes <- c("MGW", "ProT", "HelT", "Roll")
    op <- par(mfcol = c(2, 2), mar = c(4, 4, 1, 2))
    for (i in 1:length(shapeTypes)) {
        shapes <- readShape(list.files(dirname(faFile), pattern = shapeTypes[i], full.names = TRUE))
        if(!is.null(subset))
          plotShape(shapes[subset], sequences[subset], shapeTypes[i])
        else
          plotShape(shapes, sequences, shapeTypes[i])
    }
    par(op)
}

plotShapes(faFile, subset=c(101, 102))

