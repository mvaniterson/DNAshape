library(DNAshape)
library(Biostrings)

##after reinstallation set executable
system(paste("chmod ugo+x", "/home/mvaniterson/R/x86_64-unknown-linux-gnu-library/3.2/DNAshape/v2.6_noSim/prediction.exe"))

###############################################################################
##Dickerson one sequence
faFile <- file.path(path.package("DNAshape"), "v2.6_noSim/Example/Dickerson/seq.txt")
sequence <- readDNAStringSet(faFile, format="fasta")
sequence

##do shape prediction
predictDNAShapes(faFile)

##plot single shape
shapeFiles <- list.files(path=dirname(faFile), pattern="seq")[-1]
shapeFiles
op <- par(mfcol=c(2,2))
for(i in 1:4)
  {
    shape <- readShape(file.path(dirname(faFile), shapeFiles[i]))
    shapeType <- gsub("seq.txt.", "", basename(shapeFiles[i]))
    plotShape(shape, sequence, shapeType=shapeType)
  }
par(op)

##plot all shapes at once
plotShapes(faFile)

###############################################################################
##FIS multiple sequences
faFile <- file.path(path.package("DNAshape"), "v2.6_noSim/Example/FIS/seq.txt")
sequences <- readDNAStringSet(faFile, format="fasta")
sequences

predictDNAShapes(faFile)

shapeFiles <- list.files(path=dirname(faFile), pattern="seq")[-1]
shapeFiles
i <- 3
shapes <- readShape(file.path(dirname(faFile), shapeFiles[i]))
shapeType <- gsub("seq.txt.", "", basename(shapeFiles[i]))

plotShape(shapes, sequences, shapeType=shapeType)

##plot all shapes at once
plotShapes(faFile)


###############################################################################
##nature_ex multiple sequences
faFile <- file.path(path.package("DNAshape"), "v2.6_noSim/Example/nature_ex/nature_ex.txt")
sequences <- readDNAStringSet(faFile, format="fasta")

sequences

shapeFiles <- list.files(path=dirname(faFile), pattern="nature")[-1]
shapeFiles
i <- 3
shapes <- readShape(file.path(dirname(faFile), shapeFiles[i]))
shapeType <- gsub("nature_ex.txt.", "", basename(shapeFiles[i]))
shapes

shapeType

plotShape(shapes[1], sequences[1], shapeType)




predictDNAShapes(faFile)
