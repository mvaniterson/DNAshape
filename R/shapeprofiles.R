##' read shape profiles (query table)
##'
##' read shape profiles and extends the reverse complements
##' @title read shape profiles
##' @param file file name query table
##' @param kmerLength length of kmer default 5
##' @return shape profiles for all kmers
##' @author mvaniterson
##' @export
##' @importFrom Biostrings reverseComplement DNAStringSet
read.shape.profiles <- function(file=NULL, kmerLength=5) {
    
    if(missing(file))
      file <- file.path(path.package("DNAshape"), "v2.6_noSim/QueryTable.dat")

    ##Are these estimates from the simulation?
    ##Add description of the columns
    ##Add version/date stamp

    columns <- c("pentamer", ##pentamer
                 "minor_ave", "minor_sd", "minor_num",
                 "major_ave", "major_sd", "major_num",
                 "propel_ave", "propel_sd", "propel_num",
                 "slide1_ave", "slide1_sd", "slide1_num",
                 "slide2_ave", "slide2_sd", "slide2_num",
                 "roll1_ave", "roll1_sd", "roll1_num",
                 "roll2_ave", "roll2_sd", "roll2_num",
                 "twist1_ave", "twist1_sd", "twist1_num",
                 "twist2_ave", "twist2_sd", "twist2_num")

    kmerShapes <- read.table(file,
                             nrows=(4^kmerLength)/2,
                             row.names=1,
                             col.names=columns,
                             colClasses=c("character", rep("numeric", 27)))

    ##add reverse complement to improve prediction speed
    ##swap columns of the pair properties
    kmerShapesRevComp <- kmerShapes
    kmerShapesRevComp[, grep("2_", colnames(kmerShapes))] <- kmerShapes[, grep("1_", colnames(kmerShapes))]
    kmerShapesRevComp[, grep("1_", colnames(kmerShapes))] <- kmerShapes[, grep("2_", colnames(kmerShapes))]
    rownames(kmerShapesRevComp) <- as.character(reverseComplement(DNAStringSet(rownames(kmerShapes))))
    kmerShapes <- rbind(kmerShapes, kmerShapesRevComp)
    invisible(kmerShapes)
  }

##' generate Rdata object with kmer shape profiles
##'
##' generate Rdata object with kmer shape profiles
##' @title store shape profile
##' @return filename Rdata-object
##' @author mvaniterson
##' @export
store.shape.profiles <- function() {
    shape.profiles <- read.shape.profiles(kmerLength=5)
    save(shape.profiles, file=file.path(path.package("DNAshape"), "shape.profiles.RData"))
    print(paste("Create Rdata object:", file.path(path.package("DNAshape"), "shape.profiles.RData")))
  }
