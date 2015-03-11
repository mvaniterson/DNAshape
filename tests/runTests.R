
##for unit testing
test <- function(predicted, zip)
  {
    patterns <- c(minor_ave = "MGW", propel_ave = "ProT", roll_ave = "Roll", twist_ave = "HelT")
    files <- unzip(zip)
    on.exit(file.remove(files))
    for(i in 1:length(files))
      {
        current <- as.numeric(read.table(files[i], skip=1, sep=","))
        id <- which(gsub("^.*\\.txt\\.","", files[i]) == patterns)
        if(id > 2)
          target <- as.numeric(predicted[, names(patterns)[id]])
        else
          target <- as.numeric(predicted[, names(patterns)[id]])
        target <- as.numeric(na.omit(target))
        current <- as.numeric(na.omit(current))
        print(all.equal(target, current))
        print(paste("Predicted           :", paste(target, collapse=",")) )
        print(paste("Predicted (DNAshape):", paste(current, collapse=",")))
      }
  }

library(DNAshape)

data(shape.profiles)
head(shape.profiles)
tail(shape.profiles)
dim(shape.profiles)

sequence <- "TTTTTAGC"
sequence

shapes <- seq2shape(sequence, shape.profiles)
shapes

zip <- dir(pattern=paste0(sequence, ".zip"))
test(shapes, zip)

sequence <- as.character(reverseComplement(DNAString(sequence)))
sequence

shapes <- seq2shape(sequence, shape.profiles)
zip <- dir(pattern=paste0(sequence, ".zip"))
test(shapes, zip)
