################################################################################################## 
#   Purpose: Convert a vcf file to dataframe with genotypes (0,1,2)
#   Author: Tom Poorten, tom.poorten@gmail.com
################################################################################################## 
#   Details:
#     - the genotypes dataframe can be used in downstream analyses
#     - this script is meant to be run in an interactive R session
#     - genotype format:
#         0 represents homozygous reference allele
#         1 represents heterozygous
#         2 represents homozygous alternative allele

## set working directory to the path where vcf file is located
setwd("~/Dropbox/bioinfo/snps/")
rm(list=ls())

## create modified vcf file for read in
# replace path/to/vcfFile.vcf with the vcf filename
system("sed 's/^#CHROM/CHROM/' path/to/vcfFile.vcf > temp.vcfFile.reformat.vcf")
## read in the modified vcf file
vcf = read.table(file="temp.vcfFile.reformat.vcf",row.names=NULL,sep="\t",stringsAsFactors=F, header=T)
head(vcf)

## put CHROM and POS columns in object
chrom.pos = vcf[,1:2]
## create matrix of sample data only
vcf.mat = as.matrix(vcf[,c((which(colnames(vcf)=="FORMAT")+1):ncol(vcf))])
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

## write the table with genotypes to a new file - replace the filename as desired
tab = cbind(chrom.pos, vcf.mat)
write.table(tab,file="snps_5strains.select.GT.tab",sep="\t",quote=F,row.names=F)
# check that this file can be read in correctly
tab2 = read.table("snps_5strains.select.GT.tab", sep="\t", header=T, stringsAsFactors=F)
# clean up the temp file
system("rm temp.vcfFile.reformat.vcf")


#