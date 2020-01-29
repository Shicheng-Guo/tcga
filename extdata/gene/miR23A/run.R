# miR23A, miR27A, miR24-2
source("https://raw.githubusercontent.com/Shicheng-Guo/tcga/master/extdata/gene/miR23A/HeatMap.R")
x1<-rownames(miRNA)[grep("-23a",rownames(miRNA))]
x2<-rownames(miRNA)[grep("-27a",rownames(miRNA))]
x3<-rownames(miRNA)[grep("-24-2",rownames(miRNA))]
matrix<-t(miRNA[match(c(x1,x2,x3),rownames(miRNA)),])
cor(matrix)
pdf("hsa-miR23A-h3.pdf")
heatmap.3(cor(matrix))
dev.off()
pdf("hsa-miR23A-Heat.pdf")
HeatMap(cor(matrix))
dev.off()


tfbsrnaseq<-log(tfbsrnaseq+1,2)
y<-log(tfbsmiRNAseq[match(c('hsa-mir-32',"hsa-mir-27a","hsa-mir-24-2"),rownames(tfbsmiRNAseq)),]+1,2)
matrix<-data.frame(t(y),t(tfbsrnaseq),check.names = F)

pdf("hsa-miR-mRNA-FullHeat.pdf")
HeatMap(cor(matrix))
dev.off()
write.csv(cor(matrix),file="miR23A.correlation.csv",quote=F)
