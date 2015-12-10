##' plot minor groove width as function of the sequence
##'
##' Optionally one other sequence can be added
##' @title plot Minor groove width
##' @param shape predicted shapes
##' @param seq DNAstring containing sequence
##' @param add TRUE/FALSE
##' @param ... optional plotting parameters
##' @return called for its side-effect
##' @author mvaniterson
##' @export
plotMGW <- function(shape, seq, add=FALSE, ...) {    
    y <- shape$minor_ave
    y <- c(NA, NA, y, NA)
    x <- 1:length(y)
    args <- list(...)
    col <- ifelse(any(names(args) == "col"), args[["col"]], 1)    
    if(!add) {
        op <- par(mar=c(4,5,4,2))
        plot(x, y, type="b", ylab=expression("Minor groove width["*ring(A)*"]"), xaxt="n", bty="n", xlab="", ...)
        mtext(unlist(strsplit(as.character(seq), "")), side=1, at=1:length(x), line=0, col=col)
        abline(h=axTicks(2), col="grey", lty=3)
        par(op)
    }
    else {
        args <- list(...)
        col <- ifelse(any(names(args) == "col"), args[["col"]], 1)
        points(x, y, type="b", ...)
        mtext(unlist(strsplit(as.character(seq), "")), side=1, at=1:length(x), line=2,  col=col)
    }    
}


plotROLL <- function(){
}

plotPROT <- function(){
    
}

plotHELT <- function(){
    
}


    
