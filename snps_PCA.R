################################################################################################## 
#   Purpose: Perform PCA on snp dataset
#   Author: Tom Poorten, tom.poorten@gmail.com
################################################################################################## 
#   Details:
#     - input: GT table where genotypes are coded 0,1,2, ...
#     - output: PCA plot
#     - *WARNING*: sites with missing data are EXCLUDED in this analysis
#     - prcomp() is used for PCA

rm(list=ls()) 
setwd("~/Dropbox/bioinfo/snps/SNP_data_analysis/")

## read in snp data in GT format
filein = "snps_5strains.select.short.GT.tab"
snp = read.table(file = filein, sep="\t", header=T, stringsAsFactors=F)
head(snp);tail(snp);dim(snp)

## save snp table before further manipulations
snp.save = snp

## remove sites with missing data
snp = na.omit(snp)
# take a look at the data
head(snp); dim(snp)
# output the genotype table by individual
as.list(apply(snp[,-c(1,2)], 2, function(x) table((x))) )

## prep SNP table for PCA
snp.t = t(snp[,3:ncol(snp)])
colnames(snp.t) = apply(snp,1,function(x) paste(x[1],x[2]))
snp.t.scale = scale(snp.t)
sc <- attr(snp.t.scale, "scaled:scale")
# remove the invariant loci to each other - these may be present if some strains are removed (eg JEL423, JAM81)
length(which(sc==0)) # num of invariant loci
snp.t.scale = snp.t.scale[,which(sc!=0)]

## run PCA
snp.t.pc = prcomp(snp.t.scale)

## make the PCA plot
pdf(file=paste(unlist(strsplit(filein,split=".tab"))[1],".PCAplot.pdf",sep=""),height=7,width=7,onefile=T)
par(mfrow=c(1,1))

pc1var = signif(summary(snp.t.pc)$importance[2,1], 3)
pc2var = signif(summary(snp.t.pc)$importance[2,2], 3)
plot(predict(snp.t.pc)[,1], predict(snp.t.pc)[,2], pch=" ", 
     xlab=paste("PC-1 (",pc1var*100,"%)",sep=""), ylab=paste("PC-2 (",pc2var*100,"%)",sep=""), 
     main="PCA",sub=paste(ncol(snp.t.scale),"loci"))
text(predict(snp.t.pc)[,1],predict(snp.t.pc)[,2],labels=rownames(predict(snp.t.pc)),cex=.7)

dev.off()

# 
