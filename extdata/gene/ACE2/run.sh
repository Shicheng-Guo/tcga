memo="ACE2"
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/extdata/ESCC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/colon/master/extdata/UKBB50K/CRC-UKBB-RSKAT2020.csv",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/housekeeping/master/housekeeping.txt",as.is=T,head=T,sep="\t")
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholTargetRegion50.txt",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/brca/master/extdata/brca.tcga.target.hg19.bed",as.is=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/ACE2/2019nCoV.txt",as.is=T)
head(data)
symbol<-data[,1]
source("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/ACE2/Pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/OROSmerge.R")
pancancermetadge(symbol,memo)
pancancermetaOSHR(symbol,memo)
OROSmerge(memo)


