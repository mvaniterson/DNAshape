##internal function to generate sequence/step sequence
getStepSequence <- function(sequence)
  {
    sequence <- unlist(strsplit(as.character(sequence), ""))
    repseq <- rep(sequence, each=2)
    repseq <- repseq[-c(1, length(repseq))] ##remove begin and end
    seq <- c()
    for(i in seq(1, length(repseq)-1, 2))
      seq <- c(seq, paste0(repseq[i:(i+1)], collapse=""))
    seq
  }

##plot single shape of single sequence
plotShapeSingle <- function(shape, sequence, shapeType)
  {
    shapeName <- names(shape)[1]
    if(typeof(shape) != "double")
      shape <- shape[[1]]
    plot(1:length(shape), shape, type="b", xaxt="n", xlab=shapeName , ylab=shapeType, las=2)
    if(shapeType %in% c("MGW", "ProT"))
      {
        sequence <- unlist(strsplit(as.character(sequence), ""))
        axis(1, at=1:length(shape), labels=sequence)
      }
    else
      axis(1, at=1:length(shape), labels=getStepSequence(sequence))
    grid()
  }

##plot single shape of multiple sequences of equal length
plotShapeMultiple <- function(shapes, sequences, shapeType)
  {
    shapes <- as.data.frame(shapes)
    op <- par(mar = c(ncol(shapes)+2, 4, 1, 2))
    matplot(1:nrow(shapes), shapes, type="b", xaxt="n", las=2, lty=1, pch=1, ylab=shapeType, xlab="")
    axis(1, at=1:nrow(shapes), labels=FALSE, line=0, tick=TRUE)
    for(i in 1:length(sequences))
      {
        if(shapeType %in% c("MGW", "ProT"))
          labels <- sequences[[i]]
        else
          labels <- getStepSequence(sequences[[i]])
        for(j in 1:nrow(shapes))
          mtext(labels[j], side=1, at=j, line=i-1, col=i)
      }
    grid()
    par(op)
  }


##' plot one of the four DNAshape predictions for < 10 sequences of equal length
##'
##' plot one of the four DNAshape predictions
##' MWG and ProT are predicted per bp
##' Roll and HelT are predicted per step
##' @title plot single DNAshape
##' @param shape numeric vector containing the structural features
##' @param sequence DNAString
##' @param shapeType one of MGW, ProT, HelT or Roll
##' @return generates a plot
##' @author mvaniterson
plotShape <- function(shapes, sequences, shapeType=c("MGW", "ProT", "HelT", "Roll"))
  {
    if(length(sequences) == 1) ##one sequence one shape
      plotShapeSingle(shapes, sequences, shapeType)
    else if(length(sequences) <= 10) ##multiple sequences multiple shapes
      {
        if(length(unique(width(sequences))) == 1) ## all sequences of equal length
          plotShapeMultiple(shapes, sequences, shapeType)
        else ## sequences with different length
          warning("plot sequences with different length")
      }
    else
      stop("Too many sequences to draw!")
  }

##old single sequence version
## plotShape <- function(shapes, sequences, shapeType=c("MGW", "ProT", "HelT", "Roll"))
##   {

##     shapeNames <- names(shapes)

##     shape <- shape[[1]]

##     plot(1:length(shape), shape, type="b", xaxt="n", xlab=shapeName, ylab=shapeType, las=2)

##         if(shapeType %in% c("MGW", "ProT"))
##           axis(1, at=1:length(shape), labels=sequence)
##         else
##           axis(1, at=1:length(shape), labels=getStepSequence(sequence))
##     grid()
##   }




##' plot all four DNAshape prediction in one figure
##'
##'
##' @title plot all four DNAshapes
##' @param faFile file containing sequences (fasta-format)
##' @return generates a plot
##' @author mvaniterson
plotShapes <- function(faFile)
  {
    sequences <- readDNAStringSet(faFile, format="fasta")
    shapeTypes <- c("MGW", "ProT", "HelT", "Roll")
    op <- par(mfcol=c(2, 2), mar=c(4,4,1,2))
    for(i in 1:length(shapeTypes))
      {
        shapes <- readShape(list.files(dirname(faFile), pattern=shapeTypes[i], full.names=TRUE))
        plotShape(shapes, sequences, shapeTypes[i])
      }
    par(op)
  }
