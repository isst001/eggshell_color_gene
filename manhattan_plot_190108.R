# only once:
install.packages("qqman")

# each time:
library(qqman)

#qqman demo
head(gwasResults)

inf = "manhattan_chr.txt"
gwasResults <- read.table(inf, header=TRUE, sep="\t", quote="")#in_fで指定したファイルの読み込み

#データサマリー
str(gwasResults)
head(gwasResults)
tail(gwasResults)
summary(gwasResults)
as.data.frame(table(gwasResults$CHR))


#window size
windows(
width = 10, height = 5, pointsize = 12, bg = "white")

# plot region size
par(plt = c(0.10, 0.95, 0.2, 0.9))


#cex.axis > font size, suggestiveline=F > remove line, ylim > set range, cex > change point type

manhattan(gwasResults, cex = 0.5, cex.axis = 0.8, suggestiveline = log(0.01 / 15250, 10)*-1, genomewideline = F, ylim=c(0,20))


#not used
#chr1-10andZ <- subset(gwasResults, CHR==c(1,2,3,4,5,6,7,8,9,10,32))
#manhattan(chr110Z, cex = 0.8, cex.axis = 0.8, suggestiveline = F, genomewideline = F, ylim=c(0,20)) 
#chr11-31<- subset(gwasResults, CHR==c(11:31))
#manhattan(chr1131, cex = 0.8, cex.axis = 0.8, suggestiveline = F, genomewideline = F, ylim=c(0,20))
#change colour
#LG22 <- gwasResults[gwasResults$CHR == 6,]
#snpsOfInterest <- LG22$SNP
#manhattan(gwasResults, cex = 0.8, cex.axis = 0.8, highlight = snpsOfInterest, suggestiveline = -log(0.01/15250,10), genomewideline = F, ylim=c(0,20))
#LG22
#LG22 <- subset(gwasResults, CHR== 29)
#manhattan(LG22, cex.axis = 0.6, suggestiveline = F, genomewideline = F, xlim=c(0,1.5), ylim=c(0,20))
#qqplot
#qq(gwasResults$P, main="Q-Q plot of P-values", cex = 0.5)
