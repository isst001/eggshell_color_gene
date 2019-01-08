files = list.files()

#only first time
install.packages("qtl")

library(qtl)

inf = "chr6_marker_missingdata_lessthan10per_for_lodscore180824.csv"

testmap <- read.cross("csvr", file=inf, estimate.map=FALSE) 
head(testmap)

rf <- pull.rf(testmap, what=c("rf"))
lod <- pull.rf(testmap, what=c("lod"))
write.table(rf, "result_rf_eggchr6.txt", sep="\t", append=FALSE, quote=FALSE,row.names=TRUE, col.names=TRUE)
write.table(lod, "result_lod_eggchr6.txt", sep="\t", append=FALSE, quote=FALSE,row.names=TRUE, col.names=TRUE)

##not used
## data_summary
##gt <- geno.table(testmap)
## filter
##gt2 <- gt[gt$P.value > 0.05/totmar(testmap),]
##write.table(gt, "result_chisquare_chr6.txt", sep="\t", append=FALSE, quote=FALSE,row.names=TRUE, col.names=TRUE)
##todrop <- rownames(gt[gt$P.value < 0.05/totmar(testmap),])
##todrop <- todrop[-16]
##testmap2 <- drop.markers(testmap, todrop)
##rf2 <- pull.rf(testmap2, what=c("rf"))
##lod2 <- pull.rf(testmap2, what=c("lod"))
##write.table(rf2, "result_rf_eggchr6_2.txt", sep="\t", append=FALSE, quote=FALSE,row.names=TRUE, col.names=TRUE)
##write.table(lod2, "result_lod_eggchr6_2.txt", sep="\t", append=FALSE, quote=FALSE,row.names=TRUE, col.names=TRUE)
