## analyze GTEx artery data to look for variability between samples
setwd("~/Dropbox/projects/gtex/artery")
load(file="gtex-gene-counts-artery.rda")

## transform data to a DESeq object
library(DESeq2)

## split based on coronary or tibial artery
iTibial <- which(arterytab$SMTSD == "Artery - Tibial")
coronaryArteryDES <- DESeqDataSetFromMatrix(countData = arterydat[,-iTibial],
                                            colData = arterytab[-iTibial,],
                                            design = ~ 1)
rownames(coronaryArteryDES) <- gtab$Name

tibialArteryDES <- DESeqDataSetFromMatrix(countData = arterydat[,iTibial],
                                          colData = arterytab[iTibial,],
                                          design = ~ 1)
rownames(tibialArteryDES) <- gtab$Name

## variance stabilizing transformation
coronaryArteryVSD <- varianceStabilizingTransformation(coronaryArteryDES)
save(coronaryArteryVSD, gtab, file="coronary-artery-vsd-unfiltered.rda")

tibialArteryVSD <- varianceStabilizingTransformation(tibialArteryDES)
save(tibialArteryVSD, gtab, file="tibial-artery-vsd-unfiltered.rda")

##############################################################

## filter low expressed genes - coronary
m <- rowMeans(assay(coronaryArteryVSD))
pdf("coronary-avg-vsd-histogram.pdf", width=8, height=8)
hist(m, breaks=25, main="", xlab="Normalized transformed counts")
abline(v=5, lwd=3, lty=2)
dev.off()

## select clearly expressed genes
ind <- which(m>5)
gtabMeanFiltered <- gtab[ind,]
coronaryArteryVSDMeanFiltered <- coronaryArteryVSD[ind,]
save(gtabMeanFiltered, coronaryArteryVSDMeanFiltered, file="coronary-artery-vsd-mean-filtered.rda")

## filter raw count DESeq object
coronaryArteryDESMeanFiltered <- coronaryArteryDES[ind,]
identical(rownames(coronaryArteryDESMeanFiltered), gtabMeanFiltered$Name)
save(coronaryArteryDESMeanFiltered, gtabMeanFiltered, file="coronary-artery-counts-mean-filtered.rda")

##############################################################

## filter low expressed genes - tibial
m <- rowMeans(assay(tibialArteryVSD))
pdf("tibial-avg-vsd-histogram.pdf", width=8, height=8)
hist(m, breaks=25, main="", xlab="Normalized transformed counts")
abline(v=5, lwd=3, lty=2)
dev.off()

## select clearly expressed genes
ind <- which(m>5)
gtabMeanFiltered <- gtab[ind,]
tibialArteryVSDMeanFiltered <- tibialArteryVSD[ind,]
save(gtabMeanFiltered, tibialArteryVSDMeanFiltered, file="tibial-artery-vsd-mean-filtered.rda")

## filter raw count DESeq object
tibialArteryDESMeanFiltered <- tibialArteryDES[ind,]
identical(rownames(tibialArteryDESMeanFiltered), gtabMeanFiltered$Name)
save(tibialArteryDESMeanFiltered, gtabMeanFiltered, file="tibial-artery-counts-mean-filtered.rda")


