
inf <- "WGS_genotype_dataset_for_calcuration_of_hetero_rat_180108.txt"
data <- read.table(inf, sep="\t", quote="",header=T)#in_fで指定したファイルの読み込み

#position and genotype
df <- data[,c(2,9:15)]

#dataframe for output
out <- data.frame(matrix(0,1000,8))

#calculation

for (k in 1:7) {
j = 1 + k
DF <- df[c(1,j)]

for (i in 0:999) {
h <- i +1
a <- 14000001
b <- 14001000
c <- a + 1000*i
d <- b + 1000*i
e <- DF[DF[1] >=c & DF[1]< d,]
head(e)
A <- nrow(e[e[,2]=="1/1" | e[,2]=="2/2" | e[,2]=="3/3",])
B <- nrow(e[e[,2]=="0/1" | e[,2]=="0/2" | e[,2]=="0/3" | e[,2]=="1/2" | e[,2]=="1/3" | e[,2]=="2/1"| e[,2]=="2/3"| e[,2]=="3/1"| e[,2]=="3/2",])
C <- B/(A+B)*100
out[h,1] <- c
out[h,j] <- C
}
}

colnames(out) <- colnames(df)

write.table(na.omit(out), "hetrate_results_causative1000_190108.txt", sep="\t", append=FALSE, quote=FALSE,row.names=FALSE)
