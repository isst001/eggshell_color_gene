inf <- "rowcount_for_DEGanalysis_Uterus_WT_vs_Mut_190108.txt"
x <- read.table(inf, header=T, sep="\t", row.names=1, quote="")#in_fで指定したファイルの読み込み

x2 <- x[!duplicated(x$gene_short_name),]

out <- "rowcount_for_DEGanalysis_Uterus_WT_vs_Mut_unique_190108.txt"
write.table(x2,out,append=FALSE,quote=FALSE,row.names=T,col.names=T,sep="\t")

y <- round(x2[,3:14],0)
head(y)

y2 <- y[rowSums(y) > 0,]
head(y2)
dim(y2)

group <- factor(c("A","A","A","A","A","A",##WT_Ut
"B","B","B","B","B","B" ##Mut_Ut
)) 

design <- model.matrix(~ group)
design

count <- y2

library("edgeR")

d <- DGEList(counts = count, group = group)
d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

fit <- glmFit(d, design)
param_fdr <- 0.05
foldchange <- log(2,2)

#AvsB (wild-type uterus vs. mutant uterus)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
table1 <- as.data.frame(topTags(lrt, n = nrow(count)))
table1$logFCrev <- table1$logFC
head(table1)

estimatedDEG <- table1$FDR
FDR <- estimatedDEG
logFC <- table1$logFC
estimatedDEG[FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr | abs(logFC) <= foldchange] <- 0 
table1 <- cbind(table1, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & table1$logFCrev > 0] <-1
direction[estimatedDEG==1 & table1$logFCrev < 0] <-3
direction[estimatedDEG==0] <-2
table1 <- cbind(table1, direction)
head(table1)
write.table(table1, file = "DEGanalysis_result_wt_vs_mut_uterus.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)

a <- "rowcount_for_DEGanalysis_Uterus_WT_vs_Mut_unique_190108.txt"
b <- "DEGanalysis_result_wt_vs_mut_uterus.txt"

x <- read.table(a, header=T, sep="\t", quote="")#in_fで指定したファイルの読み込み
y <- read.table(b, header=T, sep="\t", quote="")#in_fで指定したファイルの読み込み
head(x)
head(y)

x$X <- rownames(x)
y$X <- rownames(y)

temp = merge (x, y, by.x = "X", by.y = "X", sort =T, all = T)

write.table(temp, file = "DEGanalysis_result_wt_vs_mut_uterus_plus_info.txt", col.names = T, row.names=FALSE,sep = "\t", quote=FALSE)
