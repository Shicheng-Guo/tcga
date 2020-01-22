##############################################################################################################
# The script was used to do co-expression analysis between miR32A and 104 TFs within miR32A promoter regions
# Shihcheng.Guo@Gmail.com
# 2020/01/22
##############################################################################################################

id2phen4<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
}

id2phen3<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*"))
}

id2bin<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
  as.numeric(lapply(strsplit(filename,"-"),function(x) x[4]))
}

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}

RawZeroRemove<-function(data,missratio=0.5){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) sum(x==0)>=threshold))
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}


load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")

phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")

phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
head(phen1)
head(phen2)
head(phen)

colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]

phen1<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen2<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen3<-id2bin(phen$cases.0.samples.0.submitter_id)

phen$bin=phen3

include<-which(c(phen$bin==1 | phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]
colnames(input)<-phen$id
Seq<-paste(phen$project_id,phen$bin,sep="-")


save(rnaseqdata,file="~/hpc/rnaseqdata.pancancer.RData")

load("~/hpc/project/TCGA/pancancer/miRNA/data/TCGA-Pancancer.miRNAseq.RData")

rnaseqdata<-input

miRNA[1:5,1:4]

ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}

symbol<-as.character(ENSG2Symbol(rownames(rnaseqdata)))

tfbs<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/miR23A/miR23A-TF.txt",head=T,sep="\t",as.is=T)

tfbsrnaseq<-rnaseqdata[unique(na.omit(match(tfbs[,1],symbol))),]
tfbsmiRNAseq<-miRNA[,na.omit(match(colnames(tfbsrnaseq),colnames(miRNA)))]
tfbsrnaseq<-tfbsrnaseq[,match(colnames(tfbsmiRNAseq),colnames(tfbsrnaseq))]
rownames(tfbsrnaseq)=as.character(ENSG2Symbol(rownames(tfbsrnaseq)))
tfbsrnaseq[1:5,1:4]
tfbsmiRNAseq[1:5,1:4]

tfbsrnaseq<-log(tfbsrnaseq+1,2)
y<-log(tfbsmiRNAseq[match('hsa-mir-32',rownames(tfbsmiRNAseq)),]+1,2)

fit<-c()
for(i in 1:nrow(tfbsrnaseq)){
  temp<-summary(lm(y~tfbsrnaseq[i,]))
  fit<-rbind(fit,temp$coefficients[2,])
}
rownames(fit)<-rownames(tfbsrnaseq)
fit<-fit[order(fit[,4]),]
write.csv(fit,file="miR32A-TFs.coexpression.csv",quote=F)

sel<-match('hsa-mir-32',rownames(tfbsmiRNAseq))
tfbsmiRNAseqmiR32A<-tfbsmiRNAseq[(sel-2):(sel+2),]
write.csv(tfbsmiRNAseqmiR32A,file="miRNA.csv",quote=F)
write.csv(tfbsrnaseq,file="mRNAseq.csv",quote=F)

x<-tfbsrnaseq[match("RUNX1T1",rownames(tfbsrnaseq)),]

pdf("RUNX1T1.pdf")
plot(y~x,pch=16,cex=0.3,col="red",xlab="RUNX1T1",ylab="miR23A")
dev.off()

x<-tfbsrnaseq[match("GTF2B",rownames(tfbsrnaseq)),]

pdf("GTF2B.pdf")
plot(y~x,pch=16,cex=0.3,col="red",xlab="GTF2B",ylab="miR23A")
dev.off()




