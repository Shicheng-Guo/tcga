setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/temp")

temp1<-read.table("Pancancer.radiationresponse.txt",sep="\t",head=T,as.is=T)
temp2<-read.table("pancancer.chemotherapy.response.txt",sep="\t",head=T,as.is=T)

head(temp1)
head(temp2)

temp1<-subset(temp1,measure_of_response !="")
temp2<-subset(temp2,measure_of_response !="")

head(temp1)
head(temp2)


temp1<-temp1[match(unique(temp1[,2]),temp1[,2]),]
temp2<-temp2[match(unique(temp2[,1]),temp2[,1]),]

head(temp1)
head(temp2)

dim(temp1)
dim(temp2)

out1<-temp1[temp1[,2] %in% temp2[,1],]
out2<-temp2[temp2[,1] %in% temp1[,2],]

dim(out1)
dim(out2)

write.csv(temp1,file="Pancancer.radiationresponse.uni.txt",quote=F)
write.csv(temp2,file="pancancer.chemotherapy.response.uni.txt",quote=F)
write.csv(out1,file="Pancancer.radiationresponse.chemotherapy.uni.csv",quote=F)
write.csv(out2,file="Pancancer.chemotherapy.radiationresponse.uni.csv",quote=F)

