##' read DNAshape prediction
##'
##' read DNAshape prediction
##' @title read DNAshape prediction
##' @param file file containing predicted shape characteristics
##' @return list with DNAshape prediction each element is the predicted shape of a input sequence
##' @author mvaniterson
readShape <- function(file)
  {
    f <- file(file, "r")
    on.exit(close(f))
    shapes <- list()
    nseqs <- 0
    while(length(line <- readLines(f, n=1)) > 0)
      {
        nseqs <- nseqs + 1
        if(grepl(">", line))
          if(nchar(line) <= 1)
            header <- nseqs
          else
            header <- line
        else
          {
          shapes[[header]] <- eval(parse(text=paste("c(", line, ")")))
        }
      }    
    shapes
  }


