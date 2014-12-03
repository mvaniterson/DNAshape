
library(BSgenome.Hsapiens.UCSC.hg19)
library(gwascat)
library(GenomicRanges)
library(DNAshape)

data(gwrngs19)
gwrngs19

catalog <- as.data.frame(gwrngs19)
intergenic <- catalog[catalog$Context=="Intergenic",]

window <- 5
gr <- GRanges(seqnames = Rle(intergenic$seqnames), ranges = IRanges(start =intergenic$start  - window, end = intergenic$end + window))
gr

risk <- intergenic$Strongest.SNP.Risk.Allele
risk <- gsub("^.*-", "", risk)

seqs <- getSeq(Hsapiens, gr)
names(seqs) <- intergenic$SNPs

ref <- as.character(subseq(seqs, start=6, end=6))

table(ref == risk)
##FALSE  TRUE 
## 4881  1678 

dataDr <- "/media/Storage/Veni/dnashape"
writeXStringSet(seqs, file.path(dataDr, "gwascat_nc.fa"))

faFile <- file.path(dataDr, "gwascat_nc.fa")
predictDNAShapes(faFile)

shapeFiles <- list.files(path=dirname(faFile), pattern="gwascat_nc")[-1]

i <- 1
shape <- readShape(file.path(dirname(faFile), shapeFiles[i]))
shapeType <- gsub("gwascat_nc.", "", basename(shapeFiles[i]))
plotShape(shape, seqs, shapeType=shapeType)



shapes <- do.call(rbind.data.frame, shape)
colnames(shapes) <- 1:10
head(shapes)

matplot(t(shapes), type="l")
