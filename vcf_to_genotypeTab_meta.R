################################################################################################## 
#   Purpose: Convert a vcf file to dataframe with genotypes (0,1,2)
#   Author: Tom Poorten, tom.poorten@gmail.com
################################################################################################## 
#   Details:
#     - input: vcf file (e.g. GATK vcf file)
#     - output: tab-delimited table with genotypes
#     - the genotypes dataframe can be used in downstream analyses
#     - this script is meant to be run in an interactive R session
#     - genotype format:
#         0 represents homozygous reference allele
#         1 represents heterozygous
#         2 represents homozygous alternative allele
#  	      NA represents missing data

## set working directory to the path where vcf file is located
setwd("~/Dropbox/bioinfo/snps/SNP_data_analysis/")
rm(list=ls())

## create modified vcf file for read in
# replace path/to/vcfFile.vcf with the vcf filename
#*****system("sed 's/^#CHROM/CHROM/' path/to/vcfFile.vcf > temp.vcfFile.reformat.vcf")
system("sed 's/^#CHROM/CHROM/' ../snps_5strains.select.short.vcf > temp.vcfFile.reformat.vcf")
## read in the modified vcf file
vcf = read.table(file="temp.vcfFile.reformat.vcf",row.names=NULL,sep="\t",stringsAsFactors=F, header=T)
head(vcf)

## put CHROM and POS columns in object
chrom.pos = vcf[,1:2]
## create matrix of sample data only
vcf.mat = as.matrix(vcf[,c((which(colnames(vcf)=="FORMAT")+1):ncol(vcf))])
# save for metadata extraction
vcf.mat.meta = vcf.mat
head(vcf.mat)

## Convert for .vcf to .tab
for(col in 1:ncol(vcf.mat)){
  print(col)
  vcf.mat[grep("0/0",vcf.mat[,col]),col] = "0"
  vcf.mat[grep("0/1",vcf.mat[,col]),col] = "1"
  vcf.mat[grep("1/1",vcf.mat[,col]),col] = "2"
  vcf.mat[grep("./.",vcf.mat[,col]),col] = NA
}
head(vcf.mat)

## extract genotype metadata - FORMAT column contains GT:AD:DP:GQ:PL
#
allele.depth = read.depth = geno.qual = phred.likel = as.data.frame(vcf.mat.meta)
for (col in 1:ncol(vcf.mat.meta)) {
  print(col)
  allele.depth[,col] = sapply(vcf.mat.meta[,col], function(x) unlist(strsplit(x,split=":"))[2])
  read.depth[,col] = as.numeric(sapply(vcf.mat.meta[,col], function(x) unlist(strsplit(x,split=":"))[3]))
  geno.qual[,col] = as.numeric(sapply(vcf.mat.meta[,col], function(x) unlist(strsplit(x,split=":"))[4]))
  phred.likel[,col] = sapply(vcf.mat.meta[,col], function(x) unlist(strsplit(x,split=":"))[5])
}
allele.depth = cbind(vcf[,1:2], allele.depth)
read.depth = cbind(vcf[,1:2], read.depth)
geno.qual = cbind(vcf[,1:2], geno.qual)
phred.likel = cbind(vcf[,1:2], phred.likel)

## write the table with genotypes to a new file - replace the filename as desired
tab = cbind(chrom.pos, vcf.mat)
write.table(tab,file="snps_5strains.select.GT.tab",sep="\t",quote=F,row.names=F)
write.table(allele.depth,file="snps_5strains.select.AD.tab",sep="\t",quote=F,row.names=F)
write.table(read.depth,file="snps_5strains.select.DP.tab",sep="\t",quote=F,row.names=F)
write.table(geno.qual,file="snps_5strains.select.GQ.tab",sep="\t",quote=F,row.names=F)
write.table(phred.likel,file="snps_5strains.select.PL.tab",sep="\t",quote=F,row.names=F)
# check that this file can be read in correctly
tab2 = read.table("snps_5strains.select.GT.tab", sep="\t", header=T, stringsAsFactors=F)
# clean up the temp file
system("rm temp.vcfFile.reformat.vcf")

#############################################################################################
#  Get summary data for metadata - eg read depth

apply(read.depth[,-c(1:2)], 2, function(x) summary(x))
apply(geno.qual[,-c(1:2)], 2, function(x) summary(x))

boxplot(read.depth[,-c(1:2)])
