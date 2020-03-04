setwd("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp")
#install.packages(c("SNPRelate", "ggplot2", "devtools", "dplyr", "reshape2"))
library(SNPRelate)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)

#replace "DPSCF3" "" -- good_pos.vcf

vcf.fn<-"good_pos.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")
snpgdsSummary("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.table("pop.group", header = F)




ccm_pca<-snpgdsPCA(genofile)
#names(ccm_pca)
head(cbind(sample.id, pop_code))

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
names(snpset)
head(snpset$chr1)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=1)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

sample.id = pca$sample.id
pop = factor(pop_code$V1)

tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code$V1)[match(pca$sample.id, sample.id)],colour=factor(pop_code$V2)[match(pca$sample.id, sample.id)] , EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
write.csv(tab, file = "PCA_final_data.csv")
tab <- read.csv(file="PCA_final_data.csv", header=TRUE, sep=",")

head(tab)

pdf("rplot.pdf") 
plot(tab$EV2~tab$EV1, col=as.character(tab$colour), xlab="Eigenvector 2", ylab="Eigenvector 1",pch=19, cex.lab=1.5)
#with(tab, text(tab$EV2~tab$EV1, labels = as.character(tab$sample.id)), pos = 4)
dev.off() 

pdf("rplot_2.pdf") 
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.character(tab$colour), labels=lbls, pch=19)
dev.off() 





