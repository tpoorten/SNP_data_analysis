################################################################################################## 
#   Purpose: Identify Loss of Heterozygosity Regions with a Hidden Markov Model analysis
#   Author: Tom Poorten, tom.poorten@gmail.com
################################################################################################## 
#   Details:
#     - input: tab-delimited table with genotypes
#     - output: LOH calls in text file and in an RData file
#     - this script is meant to be run in an interactive R session
#     - genotype format:
#         0 represents homozygous reference allele
#         1 represents heterozygous
#         2 represents homozygous alternative allele
#          NA represents missing data
#     - WARNING - the R package "RHmm" was removed from CRAN, so an archived package 
#                 needs to be used with an older version of R (< 3.0).
#                 Here I'm using R version 2.14.2 (2012-02-29)

## load HMM package
# download archived version: http://cran.us.r-project.org/src/contrib/Archive/RHmm/RHmm_2.0.1.tar.gz
# then install from the source
install.packages("~/Software/RHmm_2.0.1.tar.gz", type="source", repos=NULL) 
library(RHmm)
  
## set working directory to the path where genotype file is located
setwd("~/Dropbox/bioinfo/snps/SNP_data_analysis/")
rm(list=ls())

## read in genotype file
snp=read.table(file="snps_5strains.select.GT.tab",sep="\t",header=T, stringsAsFactors=F)
head(snp)
## run analysis on a subset for code testing
#snp = snp[,1:4]  

## Run HMMFit function for each chromosome in each stain
# set up some objects for the loop
dir.create("HMM_LOH")
chrom = sapply(snp$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chrom = as.numeric(chrom)
head(chrom)

# set window size (bp)
window=100
i = 11; j = 3 # initialize for code testing
# run the for loop through all 15 chromosomes, and write out results to files
for(i in 1:15){
  save.results = NULL
  tmp = snp[which(chrom==i),]
  # loop through all strains
  for(j in 3:ncol(snp)){
    print(j)
    colnames(snp)[j]
    # get positions of het sites
    pos = tmp$POS[which(tmp[,j]==1)]
    h.cont = sapply(1:round(tmp$POS[length(tmp$POS)]/window), function(x) length(which( (pos > (x*window - (window-1)) )  & (pos < (x*window)) ) ))       
    h.fit = HMMFit(obs = h.cont, dis="DISCRETE", nStates=2, control=list(nInit = 10))
    h.states = viterbi(h.fit, h.cont)
    h.res = h.states$states
    save.results = rbind(save.results, h.res)
  }
  write.table(save.results, paste("HMM_LOH/discrete_states_window-",window,"_contig-",i,".txt", sep=""), sep="\t", quote=F,col.names=F, row.names=F)
}

################################################################################
################################################################################
######## Plot results of the HMM 
library(IRanges)

## set working directory to the path where genotype file is located
setwd("~/Dropbox/bioinfo/snps/SNP_data_analysis/")
rm(list=ls())

## read in genotype file
snp=read.table(file="snps_5strains.select.GT.tab",sep="\t",header=T, stringsAsFactors=F)
head(snp)
## run analysis on a subset for code testing
#snp = snp[,1:4]  

chrom = sapply(snp$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chrom = as.numeric(chrom)
head(chrom)


## make function to filter out short LOH calls
filter.hmm = function(hmm.res){
  # take out short (< 100 bp) regions of het states
  hmm.res.sh = whichAsIRanges(hmm.res==1)
  hmm.res.sh = hmm.res.sh[which( width(hmm.res.sh) < 100 )]
  hmm.res2 = hmm.res
  hmm.res2[unlist(hmm.res.sh)] = 2
  # get only long LOH (> 1000 bp) regions
  hmm.res.loh = whichAsIRanges(hmm.res2==2)
  hmm.res.loh = hmm.res.loh[which( width(hmm.res.loh) > 1000 )]
  hmm.res.loh.unlist = unlist(hmm.res.loh)
  return(hmm.res.loh.unlist)
}

## Apply filter and plot the results
# set window size (bp)
window=100
i = 14; j = 3 # initialize for code testing
dir.create("HMM_LOH/plot_results/")
for(i in 1:15){
  save.results = as.matrix(read.delim(file=paste("HMM_LOH/discrete_states_window-",window,"_contig-",i,".txt", sep=""), header=F))
  filtered = apply(save.results, 1, function(x) filter.hmm(x))
  
  png(filename=paste("HMM_LOH/plot_results/discrete_filter1_window-",window,"_contig",i,".png",sep=""),width=2000,height=1000)
  tmp = snp[which(chrom==i),]
  plot(1:tmp$POS[nrow(tmp)], pch="", ylim=c(5,(ncol(tmp)*2)), ylab="Het", xlab="position",yaxt="n")
  
  for(j in 3:ncol(snp)){
    # plot the HET sites
    pos = tmp$POS[which(tmp[,j]==1)]
    points(pos, rep.int(j*2-1,times=length(pos)) , col=rgb(0,0,0,50,maxColorValue=255) , pch=15, cex=1) 
    # plot the missing data
    pos.na = tmp$POS[which(is.na(tmp[,j]))]
    points(pos.na, rep.int(j*2-0.5,times=length(pos.na)) , col=rgb(200,0,0,50,maxColorValue=255) , pch=15, cex=1) 
    # plot the LOH calls
    if(is.list(filtered)){
      points(filtered[[j-2]]*window-window/2, rep.int(j*2, times=length(filtered[[j-2]])), col="blue", cex=1, pch=15) 
    }
    # add strain name
    text(0, j*2-1, colnames(tmp)[j], cex=1, pos=2)
  }
  dev.off()
}


######################################
######################################
## Get start/stop sites of the LOHs and save R object for each chromosome
# this may not be necessary, but is useful for downstream processing 
#     - e.g.  manually clean LOH calls, make snp file with LOH calls

## make function to filter out short LOH calls
filter.hmm.start.stop = function(hmm.res){
  # take out short (< 8,000 bp) regions of het states
  hmm.res.sh = whichAsIRanges(hmm.res==1)
  hmm.res.sh = hmm.res.sh[which( width(hmm.res.sh) < 100 )]
  hmm.res2 = hmm.res
  hmm.res2[unlist(hmm.res.sh)] = 2
  # get only long LOH (> 1000 bp) regions
  hmm.res.loh = whichAsIRanges(hmm.res2==2)
  hmm.res.loh = hmm.res.loh[which( width(hmm.res.loh) > 1000 )]
  hmm.res.loh = IRanges(hmm.res.loh)
  # multiply by window
  end(hmm.res.loh) = end(hmm.res.loh)*window
  start(hmm.res.loh) = start(hmm.res.loh)*window
  return(hmm.res.loh)
}

i=15
window = 100
for(i in 1:15){
  save.results = as.matrix(read.delim(file=paste("HMM_LOH/discrete_states_window-",window,"_contig-",i,".txt", sep=""), header=F))
  nam = paste("contig",i,sep=".")
  xx = apply(save.results, 1, filter.hmm.start.stop)
  names(xx) = colnames(snp)[3:ncol(snp)]
  assign(nam, xx)
}
save(contig.1, contig.2, contig.3, contig.4, contig.5, contig.6, contig.7, contig.8, contig.9, contig.10, contig.11, contig.12, contig.13, contig.14, contig.15, file="HMM_LOH/LOH_start_stops_window-100.RData")

# check to make sure that the R object saved correctly
rm(list=ls())
library(IRanges)
load(file="HMM_LOH/LOH_start_stops_window-100.RData")
# check out contig 1
contig.1
#####################################################################################################
#####################################################################################################
# Clean up the breakpoints manually if needed (eg. if a known centromere is present)

rm(list=ls())
library(IRanges)
load(file="HMM_LOH/LOH_start_stops_window-100.RData")

########
# Contig 1
#
ii = 2
contig.1[[ii]] = contig.1[[ii]][-c(1)]

# contig 2
ii = 2
contig.2[[ii]] = contig.2[[1]]

# save the cleaned HMM calls
save(contig.1, contig.2, contig.3, contig.4, contig.5, contig.6, contig.7, contig.8, contig.9, contig.10, contig.11, contig.12, contig.13, contig.14, contig.15, file="HMM_LOH/LOH_start_stops_window-100_CLEANED.RData")
##############################################
##############################################
# Plot the cleaned breakpoints
library(IRanges)
rm(list=ls())
#load(file="HMM_LOH/LOH_start_stops_window-100.RData")
load(file="HMM_LOH/LOH_start_stops_window-100_CLEANED.RData")


## read in genotype file
snp = read.table(file="snps_5strains.select.GT.tab",sep="\t",header=T, stringsAsFactors=F)

## run analysis on a subset for code testing
#snp = snp[,1:4]  

## prep objects for plotting
contig.obj = ls()[grep("contig",ls())]
reorder = as.numeric(sapply(contig.obj, function(x) unlist(strsplit(x, split="\\."))[2]))
contig.obj = contig.obj[order(reorder)]
contig.obj

chrom = sapply(snp$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chrom = as.numeric(chrom)

## plot the clean LOH calls from the saved R objects
i=14
for(i in 1:15){
  print(i)
  png(filename=paste("HMM_LOH/plot_results/CLEANED_contig",i,".png",sep=""),width=2000,height=1000)
  tmp = snp[which(chrom==i),]
  plot(1:tmp$POS[nrow(tmp)], pch="", ylim=c(5,(ncol(tmp)*2)), ylab="Het", xlab="position",yaxt="n")
  for(j in 3:ncol(snp)){
    pos = tmp$POS[which(tmp[,j]==1)]
    #    points(pos, rep.int(j*2-1,times=length(pos)) , col=rgb(0,0,0,50,maxColorValue=255) , pch=15, cex=.4)
    points(pos, rep.int(j*2-1,times=length(pos)) , pch=15, cex=1, col=rgb(0,0,0,50,maxColorValue=255))
    text(0, j*2-1, colnames(tmp)[j], cex=1, pos=2)
    
    if(length(start(get(contig.obj[i])[[j-2]])) > 0){
      segments(x0 = start(get(contig.obj[i])[[j-2]]), x1 = end(get(contig.obj[i])[[j-2]]), y0 = j*2-.5 , col="blue", lwd=2)
    }
    # plot the missing data
    pos.na = tmp$POS[which(is.na(tmp[,j]))]
    points(pos.na, rep.int(j*2-0.5,times=length(pos.na)) , col=rgb(200,0,0,50,maxColorValue=255) , pch=15, cex=1) 
  }
  dev.off()
}


############################################################################
############################################################################
############################################################################
#  Remove the LOH regions
## set working directory to the path where genotype file is located
setwd("~/Dropbox/bioinfo/snps/SNP_data_analysis/")
rm(list=ls())

## read in genotype file
snp=read.table(file="snps_5strains.select.GT.tab",sep="\t",header=T, stringsAsFactors=F)
head(snp)
## run analysis on a subset for code testing
#snp = snp[,1:4]  

chrom = sapply(snp$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chrom = as.numeric(chrom)

# use only chromosomes 1 thru 15 in the snp object
snp = snp[which(chrom %in% c(1:15)),]

#load(file="HMM_LOH/LOH_start_stops_window-100.RData")
load(file="HMM_LOH/LOH_start_stops_window-100_CLEANED.RData")

# check to make sure the strains are in the right order
names(contig.2)
colnames(snp)[-c(1,2)]
# match(names(contig.2), colnames(snp)[-c(1,2)])
# snp = snp[,match(names(contig.2), colnames(snp))]

# save another snp object - for snp table with LOH genotypes recoded as "9"
snp2 = snp

# initialize for code testing
i = 2 # chrom
j = 2 # strain

## reorder the contig.obj vector
contig.obj = ls()[grep("contig",ls())]
reorder = as.numeric(sapply(contig.obj, function(x) unlist(strsplit(x, split="\\."))[2]))
contig.obj = contig.obj[order(reorder)]
contig.obj

for(i in 1:15){
  print(i)
  for (j in 1:(ncol(snp)-2)){
    snp2[unlist(sapply(1:length(get(contig.obj[i])[[j]]), function(x) which( (chrom == i) & (snp$POS > start(get(contig.obj[i])[[j]])[x] ) & (snp$POS < end(get(contig.obj[i])[[j]])[x]) ))), j+2] = 9
  }
}

write.table(snp2,file="HMM_LOH/snps_5strains.select.GT.noLOH.tab",sep="\t",quote=F,row.names=F)

# save another snp object - for snp table with non-LOH genotypes recoded as "9"
snp3 = snp
for(c in 3:ncol(snp)){
  snp3[which(snp2[,c] != 9),c] = 9
}
write.table(snp3,file="HMM_LOH/snps_5strains.select.GT.LOHonly.tab",sep="\t",quote=F,row.names=F)

