inf="raw_count_data_all_organs_tissues_190108.txt"
x <- read.table(inf, header=T, sep="\t", row.names=1, quote="")

x <- x[!duplicated(x$gene_short_name),]
write.table(x, "raw_count_data_all_organs_tissues_unique_190108.txt", row.names=T,col.names = T, sep = "\t", quote=FALSE)

inf2="raw_count_data_all_organs_tissues_unique_190108.txt"
x2 <- read.table(inf2, header=T, sep="\t", row.names=1, quote="")

y <- round(x2[,3:56],0)

y2 <- y[rowSums(y) > 0,]

group <- factor(c("A","A","A","A","A","A", ##Uterus
"B","B","B","B","B","B", ##Isthmus
"C","C","C","C","C","C", ##Brain
"D","D","D","D","D","D", ##Heart
"E","E","E","E","E","E",  ##Intestine
"F","F","F","F","F","F",  ##Kidney
"G","G","G","G","G","G",  ##Liver
"H","H","H","H","H","H", ##Lung
"I","I","I","I","I","I" ##Muscle
))

design <- model.matrix(~ group)

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

#AvsB (Uterus vs. Isthmus)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
table1 <- as.data.frame(topTags(lrt, n = nrow(count)))
table1$logFCrev <- table1$logFC*-1
head(table1)

estimatedDEG <- table1$FDR
FDR <- table1$FDR
logFC <- table1$logFCrev 
estimatedDEG[FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr | abs(logFC) <= foldchange] <- 0 
table1 <- cbind(table1, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table1 <- cbind(table1, direction)
head(table1)
write.table(table1, file = "DEGanalysis_result_Uterus_vs_Isthmus.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)

#AvsC(Uterus vs. Brain)
lrt <- glmLRT(fit, coef=3)
topTags(lrt)
table2 <- as.data.frame(topTags(lrt, n = nrow(count)))
table2$logFCrev <- table2$logFC*-1

FDR<-table2$FDR
logFC <- table2$logFCrev

estimatedDEG <- FDR
estimatedDEG[table2$FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr | abs(logFC) <= foldchange] <- 0 
table2 <- cbind(table2, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table2 <- cbind(table2, direction)
head(table2)
write.table(table2, file = "DEGanalysis_result_Uterus_vs_Brain.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)

#AvsD (Uterus vs. Heart)
lrt <- glmLRT(fit, coef = 4)
topTags(lrt)
table3 <- as.data.frame(topTags(lrt, n = nrow(count)))
table3$logFCrev <- table3$logFC*-1

FDR<-table3$FDR
logFC <- table3$logFCrev

estimatedDEG <- FDR
estimatedDEG[table3$FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr | abs(logFC) <= foldchange] <- 0 
table3 <- cbind(table3, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table3 <- cbind(table3, direction)
head(table3)

write.table(table3, file = "DEGanalysis_result_Uterus_vs_Heart.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)

#AvsE (Uterus vs. Intestine)
lrt <- glmLRT(fit, coef = 5)
topTags(lrt)
table4 <- as.data.frame(topTags(lrt, n = nrow(count)))
table4$logFCrev <- table4$logFC*-1

FDR<-table4$FDR
logFC <- table4$logFCrev

estimatedDEG <- FDR
estimatedDEG[table4$FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr | abs(logFC) <= foldchange] <- 0 
table4 <- cbind(table4, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table4 <- cbind(table4, direction)
head(table4)
write.table(table4, file = "DEGanalysis_result_Uterus_vs_Intestine.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)

#AvsF (Uterus vs. Intestine)
lrt <- glmLRT(fit, coef = 6)
topTags(lrt)
table5 <- as.data.frame(topTags(lrt, n = nrow(count)))
table5$logFCrev <- table5$logFC*-1

FDR<-table5$FDR
logFC <- table5$logFCrev

estimatedDEG <- FDR
estimatedDEG[table5$FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr| abs(logFC) <= foldchange] <- 0 
table5 <- cbind(table5, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table5 <- cbind(table5, direction)
head(table5)
write.table(table5, file = "DEGanalysis_result_Uterus_vs_Kidney.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)


#AvsG (Uterus vs. Liver)
lrt <- glmLRT(fit, coef = 7)
topTags(lrt)
table6 <- as.data.frame(topTags(lrt, n = nrow(count)))
table6$logFCrev <- table6$logFC*-1

FDR<-table6$FDR
logFC <- table6$logFCrev

estimatedDEG <- FDR
estimatedDEG[table6$FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr| abs(logFC) <= foldchange] <- 0 
table6 <- cbind(table6, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table6 <- cbind(table6, direction)
head(table6)
write.table(table6, file = "DEGanalysis_result_Uterus_vs_Liver.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)


#AvsH (Uterus vs. Lung)
lrt <- glmLRT(fit, coef = 8)
topTags(lrt)
table7 <- as.data.frame(topTags(lrt, n = nrow(count)))
table7$logFCrev <- table7$logFC*-1

FDR<-table7$FDR
logFC <- table7$logFCrev

estimatedDEG <- FDR
estimatedDEG[table7$FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr| abs(logFC) <= foldchange] <- 0 
table7 <- cbind(table7, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table7 <- cbind(table7, direction)
head(table7)
write.table(table7, file = "DEGanalysis_result_Uterus_vs_Lung.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)


#AvsI (Uterus vs. Muscle)
lrt <- glmLRT(fit, coef = 9)
topTags(lrt)
table8 <- as.data.frame(topTags(lrt, n = nrow(count)))
table8$logFCrev <- table8$logFC*-1

FDR<-table8$FDR
logFC <- table8$logFCrev

estimatedDEG <- FDR
estimatedDEG[FDR < param_fdr & abs(logFC) > foldchange] <- 1
estimatedDEG[FDR >= param_fdr | abs(logFC) <= foldchange] <- 0 
table8 <- cbind(table8, estimatedDEG)

direction <- estimatedDEG
direction[estimatedDEG==1 & logFC > 0] <-1
direction[estimatedDEG==1 & logFC < 0] <-3
direction[estimatedDEG==0] <-2
table8 <- cbind(table8, direction)
head(table8)
write.table(table8, file = "DEGanalysis_result_Uterus_vs_Muscle.txt", col.names = T, row.names = T, sep = "\t", quote=FALSE)


X <- rownames(table1)
table1 <- cbind(table1, X)

X <- rownames(table2)
table2 <- cbind(table2, X)

X <- rownames(table3)
table3 <- cbind(table3, X)

X <- rownames(table4)
table4 <- cbind(table4, X)

X <- rownames(table5)
table5 <- cbind(table5, X)

X <- rownames(table6)
table6 <- cbind(table6, X)

X <- rownames(table7)
table7 <- cbind(table7, X)

X <- rownames(table8)
table8 <- cbind(table8, X)

x<-table1 #causion
y<-table2

#Merge by X
temp1 = merge (x, y, by.x = "X", by.y = "X", sort =T, all = T)

y <-table3 #be careful
temp1 = merge (temp1, y, by.x = "X", by.y = "X", sort =T, all = T)

y<-table4 #be careful
temp1 = merge (temp1, y, by.x = "X", by.y = "X", sort =T, all = T)

y<-table5 #be careful
temp1 = merge (temp1, y, by.x = "X", by.y = "X", sort =T, all = T)

y<-table6 #be careful
temp1 = merge (temp1, y, by.x = "X", by.y = "X", sort =T, all = T)

y<-table7 #be careful
temp1 = merge (temp1, y, by.x = "X", by.y = "X", sort =T, all = T)

y<-table8 #be careful
temp1 = merge (temp1, y, by.x = "X", by.y = "X", sort =T, all = T)

#high(1), even(2), low(3) in uterus compared to other organs or tissues 
direction <- temp1[,c(9,17,25,33,41,49,57,65)]
head(direction)

concate <- paste(direction[,1], direction[,2], sep="")
concate <- paste(concate, direction[,3], sep="")
concate <- paste(concate, direction[,4], sep="")
concate <- paste(concate, direction[,5], sep="")
concate <- paste(concate, direction[,6], sep="")
concate <- paste(concate, direction[,7], sep="")
concate <- paste(concate, direction[,8], sep="")

temp2 <- cbind(temp1, concate)

write.table(temp2, file = "DEGanalysis_result_all_tissues_organs_merge.txt", col.names = T, sep = "\t", quote=FALSE)


inf <- "raw_count_data_all_organs_tissues_unique_190108.txt"
inf2 <- "DEGanalysis_result_all_tissues_organs_merge.txt"

x <- read.table(inf, header=T, sep="\t", quote="")
y <- read.table(inf2, header=T, sep="\t", quote="")

x$X <- row.names(x)

temp2 = merge (x, y, by.x = "X", by.y = "X", sort =T, all = T)

head(temp2)
write.table(temp2, file = "DEGanalysis_result_all_tissues_organs_merge_info.txt", col.names = T, row.names=FALSE,sep = "\t", quote=FALSE)

