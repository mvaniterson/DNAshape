##R-3.1.2 on VM
library(gwascat)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(DNAshape)

data(gwrngs19)
gwrngs19

catalog <- as.data.frame(gwrngs19)

intergenic <- catalog[catalog$Intergenic == "1",] ##assuming these are all noncoding
dim(intergenic)
##[1] 6559   40

head(intergenic)

intergenic$Strongest.SNP.Risk.Allele <- gsub("^.*-", "", intergenic$Strongest.SNP.Risk.Allele)
intergenic <- subset(intergenic,  Strongest.SNP.Risk.Allele %in% c("A", "C", "G", "T")) ##only single nucleotide changes
dim(intergenic)
##[1] 4727   40

head(intergenic)

window <- 5 ##gives 11 nucl. ~one pitch
gr <- GRanges(seqnames = Rle(intergenic$seqnames), ranges = IRanges(start =intergenic$start  - window, end = intergenic$end + window))
gr

intergenic$Ref.Seq <- as.character(getSeq(Hsapiens, gr))
intergenic$Ref.Allele <- substr(intergenic$Ref.Seq, 6, 6)

intergenic <- subset(intergenic, Ref.Allele != Strongest.SNP.Risk.Allele)
dim(intergenic)
##[1] 3049   42

head(intergenic)

intergenic$Alt.Seq <- intergenic$Ref.Seq
substr(intergenic$Alt.Seq, 6, 6) <- intergenic$Strongest.SNP.Risk.Allele

head(intergenic)

##predict referenc shape profiles
seqs <- DNAStringSet(intergenic$Ref.Seq)
faFile <- "ref.fa"
writeXStringSet(seqs,  faFile)
predictDNAShapes(faFile)

seqs <- DNAStringSet(intergenic$Alt.Seq)
faFile <- "alt.fa"
writeXStringSet(seqs,  faFile)
predictDNAShapes(faFile)

shapeDistances <- function(faFile1, faFile2)
  {
    shapeFiles1 <- list.files(pattern=faFile1, full.names=TRUE)[-1]
    shapeFiles2 <- list.files(pattern=faFile2, full.names=TRUE)[-1]
    d <- c()
    for(i in 1:4)
      {
        shape1 <- readShape(shapeFiles1[i])
        shape1 <- do.call(rbind.data.frame, shape1)
        colnames(shape1) <- 1:ncol(shape1)
        shape1 <- shape1[, !is.na(shape1[1,])]

        shape2 <- readShape(shapeFiles2[i])
        shape2 <- do.call(rbind.data.frame, shape2)
        colnames(shape2) <- 1:ncol(shape2)
        shape2 <- shape2[, !is.na(shape2[1,])]

        d <- cbind(d, sqrt(rowSums((shape1 - shape2)^2)))
      }
    colnames(d) <- gsub("^.*\\.", "", shapeFiles1)
    d
  }

d <- shapeDistances("ref.fa", "alt.fa")

dim(intergenic)
dim(d)

pairs(d)

intergenic <- cbind(intergenic, d)

intergenic <- intergenic[order(intergenic$HelT, intergenic$MGW, intergenic$ProT, intergenic$Roll, decreasing=TRUE),]

h <- head(intergenic, n=10)

table(h$Disease.Trait)

op <- par(mfcol=c(4, 2), mar=c(2,2,2,2))
for(diseasetrait in unique(h$Disease.Trait))
  {
    print(diseasetrait)
    s <- subset(intergenic, Disease.Trait == diseasetrait)
    print(dim(s))
    plot(s$Pvalue_mlog, main=diseasetrait, type="h", ylim=c(5, 20))
  }
par(op)


plot(intergenic$Pvalue_mlog, intergenic$OR.or.beta, log="xy")

plot(intergenic$Pvalue_mlog, intergenic$MGW, xlim=c(0, 20))
y <- sqrt(rowSums(intergenic[, 44:47]^2))
thr <- 17
text(intergenic$Pvalue_mlog[y > thr], intergenic$MGW[y > thr], intergenic$SNPs[y > thr], col=2)
intergenic$SNPs[y>thr]
