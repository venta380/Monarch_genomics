setwd("/scratch/vt20265/")
#install.packages(c("SNPRelate", "ggplot2", "devtools", "dplyr", "reshape2"))
library(SNPRelate)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)


vcf.fn<-"SNPs_PASS_only_for_1977_PCA.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")
snpgdsSummary("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.ids <- data.frame(sample.id)

pop_code <- read.table("pop.group", header = F)
df3 = merge(sample.ids, pop_code, by.x=c("sample.id"), by.y=c("V1"))



#ccm_pca<-snpgdsPCA(genofile)
#names(ccm_pca)
head(cbind(sample.id, pop_code))

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, num.thread=3)
names(snpset)
head(snpset$chr1)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

sample.id = pca$sample.id
pop = factor(pop_code$V2)

#tab <- data.frame(sample.id = ccm_pca$sample.id, pop = factor(pop_code$V1)[match(ccm_pca$sample.id, sample.id)],colour=factor(pop_code$V2)[match(ccm_pca$sample.id, sample.id)] , EV1 = ccm_pca$eigenvect[,1], EV2 = ccm_pca$eigenvect[,2], stringsAsFactors = FALSE)

tab <- data.frame(sample.id = pca$sample.id, pop = df3$V2, colour=df3$V3 , EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
write.csv(tab, file = "PCA_final_data_1977.csv")
tab <- read.csv(file="PCA_final_data_1977.csv", header=TRUE, sep=",")
#tab=tab[-c(21,22,47,48,13,17,18), ]
head(tab)

pdf("rplot2.pdf") 
plot(tab$EV2~tab$EV1, col=as.character(tab$colour), xlab="Eigenvector 2", ylab="Eigenvector 1",pch=19, cex.lab=1.5)
#with(tab, text(tab$EV2~tab$EV1, labels = as.character(tab$sample.id)), pos = 4)
dev.off() 

pdf("rplot_2.pdf") 
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.character(tab$colour), labels=lbls, pch=19)
dev.off() 



####
#LEA

library(RColorBrewer)

library(LEA)

jBrewColors <- brewer.pal(n = 10, name = "Paired")
V3<-c("02_77","03_77","04_77","05_77","07_77","08_77","10_77","12_77","13_77","14_77","ESB1","ESB10","ESB3","ESB4","ESB5","ESB8","ESB9","HH1","HH10","HH2","HH3","HH4","HH5","HH6","HH7","HH8","HH9","HI023","HI033","NJ1","NJ116","NJ203","PL1","PL10","PL2","PL4","PL5","PL6","PL7","PL8","PL9","T14","T9","mex1527","mex536","mex915","mex919","mex986","stm146","stm163")

West<-data.frame('sample.id'=c("HH1","HH2","HH3","HH4","HH5","HH6","HH7","HH8","HH9","HH10","PL1","PL2","PL4","PL5","PL6","PL7","PL8","PL9","PL10","ESB1","ESB3","ESB4","ESB5","ESB8","ESB9","ESB10"))
East<-data.frame('sample.id'=c("stm163","stm146","T9","T14","NJ203","NJ116","NJ1","HI023","HI033","mex986","mex919","mex915","mex536","mex1527"))
new_1977<-data.frame('sample.id'=c("02_77","03_77","04_77","05_77","07_77","08_77","10_77","12_77","13_77","14_77"))


#obj.snmf = snmf("new.geno", K = 1:4, ploidy = 2, entropy = T,alpha = 100, project = "new")


obj.snmf  = load.snmfProject("plink.snmfProject")

pdf(file= './error.pdf' ,onefile=T,paper='A4', )
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
dev.off() 



qmatrix1 = Q(obj.snmf, K = 2)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3= V3)

qmat2 <- qmatrix1[order(qmatrix1[,1]),]

West_2=merge(qmat2, West, by.x=c("V3"), by.y=c("sample.id"))
East_2=merge(qmat2, East, by.x=c("V3"), by.y=c("sample.id"))
new_1977_2=merge(qmat2, new_1977, by.x=c("V3"), by.y=c("sample.id"))

West_2=West_2[order(West_2[,2]),]
East_2=East_2[order(East_2[,2]),]
new_1977_2=new_1977_2[order(new_1977_2[,2]),]

West_2=data.frame(V1=West_2$V1, V2=West_2$V2)
East_2=data.frame(V1=East_2$V1, V2=East_2$V2)
new_1977_2=data.frame(V1=new_1977_2$V1, V2=new_1977_2$V2)

jBrewColors <- c('#365deb','#FF0000')
pdf(file= './K_2.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,9), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(0.3467,1.0))
plot1=barplot(t(West_2), col = jBrewColors, border = NA, ylab = "K=2", xlab = "West")
par(pin=c(0.1867,1.0))
plot1=barplot(t(East_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "East")
par(pin=c(0.0933,1.0))
plot1=barplot(t(new_1977_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "1977")
dev.off() 




qmatrix1 = Q(obj.snmf, K = 3)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3 = qmatrix1[,3], V4= V3)

qmat2 <- qmatrix_1[with(qmatrix_1, order(qmatrix_1[,3]+qmatrix_1[,2])),]
qmat2 <- qmat2[with(qmat2, order(qmat2[,1]+qmat2[,2])),]

qmat_2 = data.frame(V3 = qmat2[,3],V2 = qmat2[,2], V1 = qmat2[,1])



West_2=merge(qmat2, West, by.x=c("V4"), by.y=c("sample.id"))
East_2=merge(qmat2, East, by.x=c("V4"), by.y=c("sample.id"))
new_1977_2=merge(qmat2, new_1977, by.x=c("V4"), by.y=c("sample.id"))


West_2=West_2[with(West_2, order(West_2[,4]+West_2[,3])),]
East_2=East_2[with(East_2, order(East_2[,4]+East_2[,3])),]
new_1977_2=new_1977_2[with(new_1977_2, order(new_1977_2[,4]+new_1977_2[,3])),]

West_2=West_2[with(West_2, order(West_2[,2]+West_2[,3])),]
East_2=East_2[with(East_2, order(East_2[,2]+East_2[,3])),]
new_1977_2=new_1977_2[with(new_1977_2, order(new_1977_2[,2]+new_1977_2[,3])),]



West_2=data.frame(V1=West_2$V1, V2=West_2$V2, V3=West_2$V3)
East_2=data.frame(V1=East_2$V1, V2=East_2$V2, V3=East_2$V3)
new_1977_2=data.frame(V1=new_1977_2$V1, V2=new_1977_2$V2, V3=new_1977_2$V3)


jBrewColors <- c('#FF0000','#365deb','#7F7F7F')
pdf(file= './K_3.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,9), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(0.3467,1.0))
plot1=barplot(t(West_2), col = jBrewColors, border = NA, ylab = "K=2", xlab = "West")
par(pin=c(0.1867,1.0))
plot1=barplot(t(East_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "East")
par(pin=c(0.0933,1.0))
plot1=barplot(t(new_1977_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "1977")
dev.off() 



pdf(file= './1977_K2.pdf' ,onefile=T,paper='A4', )
plot1=barplot(t(qmat2), col = jBrewColors, border = NA, ylab = "K=3")
abline(h=0)
text(t(qmat2), labels = qmat2[,3]))
dev.off() 

qmatrix1 = Q(obj.snmf, K = 2)
pdf(file= './1977_K2.pdf' ,onefile=T,paper='A4', )
plot1=barplot(t(qmatrix1), col = jBrewColors, border = NA, ylab = "K=3")
abline(h=0)
text(t(qmat2), labels = qmat2[,3]))
dev.off() 





