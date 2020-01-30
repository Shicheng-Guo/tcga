## METTL
setwd("/home/guosa/hpc/project/esopheagul")
setwd("/home/guosa/hpc/project/crc")
setwd("/mnt/bigdata/Genetic/Projects/shg047/temp/housekeeping")
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/extdata/ESCC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/colon/master/extdata/UKBB50K/CRC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/METTL/METTL.txt",as.is=T,head=F,sep="\t")
symbol<-data[,1]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/Pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/OROSmerge.R")
memo="METTL"
pancancermetadge(symbol,memo)
pancancermetaOSHR(symbol,memo)
OROSmerge(memo)
