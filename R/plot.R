.plotShape <- function(x, y, ylab, seq, add=FALSE, ...) {
    args <- list(...)
    col <- ifelse(any(names(args) == "col"), args[["col"]], "black")
    if(!add) {
        op <- par(mar=c(4,5,4,2))
        plot(x, y, type="b", ylab=ylab, xaxt="n", bty="n", xlab="", ...)
        mtext(unlist(strsplit(as.character(seq), "")), side=1, at=1:length(x), line=0, col=col)
        abline(h=axTicks(2), col="grey", lty=3)
        par(op)
    }
    else {       
        points(x, y, type="b", ...)
        mtext(unlist(strsplit(as.character(seq), "")), side=1, at=1:length(x), line=2,  col=col)
    }
}

##' plot shape as function of function
##'
##' Optionally one other sequence can be added
##' @title plot shape as function of sequences
##' @param shape predicted shapes
##' @param seq DNAstring containing sequence
##' @param type "mgw", "prot", "roll" or "helt"
##' @param add TRUE/FALSE
##' @param ... optional plotting parameters
##' @return called for its side-effect
##' @author mvaniterson
##' @export
plotShape <- function(shape, seq, type=c("mgw", "prot", "roll", "helt"), add=FALSE, ...){
    type <- match.arg(type)
    switch(type,
           mgw = .plotShape(x = seq(length=nrow(shape)+3),
               y = c(NA, NA, shape$minor_ave, NA),
               ylab = expression("Minor groove width["*ring(A)*"]"), seq=seq, add=add, ...),
           prot = .plotShape(x = seq(length=nrow(shape)+3),
               y = c(NA, NA, shape$propel_ave, NA),
               ylab = expression("Propeller Twist["*degree*"]"), seq=seq, add=add, ...),
           roll = .plotShape(x = seq(from=0.5, length=nrow(shape)+2),
               y = c(NA, NA, shape$roll_ave),
               ylab = expression("Roll["*degree*"]"), seq=seq, add=add, ...),
           helt = .plotShape(x = seq(from=0.5, length=nrow(shape)+2),
               y = c(NA, NA, shape$twist_ave),
               ylab = expression("Helix Twist["*degree*"]"), seq=seq, add=add, ...))
}
