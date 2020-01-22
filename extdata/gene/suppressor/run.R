
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/OROSmerge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
memo="DRTSG"
data<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/suppressor/figure/TSG.pancancer.dge.tcga.pancancer.smd.meta.pvalue.csv")
int<-subset(data,beta<0 & pval<2*10^-7)
symbolist<-as.character(int[,ncol(int)])
pancancermetadge(symbolist,memo=paste(memo,".pancancer.dge",sep=""))
pancancermetaOsHr(symbolist,memo=paste(memo,".pancancer.hr.os",sep=""))
OROSmerge(memo)
