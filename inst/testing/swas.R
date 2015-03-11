library(snpStats)

chr <- 10
dataDr <- "/virdir/Scratch/mviterson/regulation/genotypes"
biobank <- "LLS"
date <- "20140402"

##get sample file
sampleFile <- paste0(biobank, "-imputed-", date, ".sample")
srm.path <- file.path("https://fly1.grid.sara.nl:2882/pnfs/grid.sara.nl/data/bbmri.nl/RP3/GWAS_ImputationGoNLv5", biobank)
SRM2VM(srm.path, sampleFile, dataDr, proxy = "/tmp/x509up_u34722")
samples <- read.table(file.path(dataDr, sampleFile), header=TRUE)

##get imputed file
snpFile <- paste0(biobank, "-imputed-chr", chr,"-", date)
SRM2VM(srm.path, snpFile, dataDr, proxy = "/tmp/x509up_u34722")
snpFile <- file.path(dataDr, snpFile)

##get info file
infoFile <- paste0(snpFile, "_info")
SRM2VM(srm.path, infoFile, dataDr, proxy = "/tmp/x509up_u34722")
infoFile <- file.path(dataDr, infoFile)
info <- read.table(infoFile)

head(info)

##remove indels
cmd <- paste("grep '[^ ]* [^ ]* [^ ]* [ACGT] [ACGT][^ACGT]'", snpFile, ">", paste(snpFile, "noindels", sep="-"))
print(system(cmd, intern=TRUE))
snpFile <- paste(snpFile, "noindels", sep="-")

##read genotypes
G <- read.impute(snpFile, snpcol=2)

g <- G[,101]
plotUncertainty(g)
g

(p <- g2post(g)) ## Transform to probabilities ...
(m <- p[,2]+2*p[,3]) ## Posterior expectations
(mg <- mean2g(m) ## Transform to codes ...)



##filter on MAF and certain.calls
MAF <- 0.05
Ccalls <- 0.4
Gsummary <- col.summary(G)

dim(G)
G <- G[,Gsummary$MAF > MAF & Gsummary$Certain.calls > Ccalls]
dim(G)


##get references sequences around snps
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

first2Rows <- read.table(snpFile, nrows=2)
ncols <- ncols(first2Rows)
##extract Annotation
cc <- rep("NULL", ncols)
cc[1:5] <- "character"
snpAnnotation <- read.table(snpFile, colClasses=cc, comment.char="")
snpAnnotation <- snpAnnotation[, -1] ##first row contains data type information
colnames(snpAnnotation) <- c("ID", "pos", "allele1", "allele2")
head(snpAnnotation)
snpAnnotation[101,]
Gsummary[101,]
 
colnames(p) <- c(paste0(snpAnnotation[101,c(3,3)], collapse=""),
              paste0(snpAnnotation[101,3:4], collapse=""),
              paste0(snpAnnotation[101,c(4,4)], collapse=""))
 
window <- 5 ##gives 11 nucl. ~one pitch
gr <- GRanges(seqnames = Rle(paste0("chr", chr)), ranges = IRanges(start = as.integer(snpAnn$pos) - window, end=as.integer(snpAnn$pos) + window))
names(gr) <- snpAnn$ID

gr
seqs <- getSeq(Hsapiens, gr)
names(seqs) <- names(gr)

##
left <- narrow(seqs, end=5)
right <- narrow(seqs, start=7)
##calculate shape profile
library(DNAshape)

for(base in c("A", "T", "C", "G"))
  {
    writeXStringSet(xscat(left, base, right), paste0("tmp", base, ".fa"))
    predictDNAShapes(paste0("tmp", base, ".fa"))
  }


g <- G[1:10, 1]

p <- g2post(g) ## Transform to probabilities ...
pg <- post2g(p) ## ... and back to codes
m <- p[,2]+2*p[,3] ## Posterior expectations
mg <- mean2g(m) ## Transform to codes ...
pmg <- g2post(mg) ## ... and transform to probabilities
## Write everything out
print(cbind(as(g, "numeric"), p, as.numeric(pg), m, as.numeric(mg), pmg))

as(g, "character")


