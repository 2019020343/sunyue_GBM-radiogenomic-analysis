#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

rm(list = ls())

library(limma)
TCGA1<-read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\HT_HG-U133A\\HT_HG-U133A.txt",header=T,sep = "\t", quote = "")
group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
less<-group_name[which(group_name$group=="LOW"),1]
more<-group_name[which(group_name$group=="HIGH"),1]
length(more)
length(less)

moredataout<-TCGA1[,1]
for(i in 1:length(more)){
 mores<-unlist(strsplit(as.character(more[i]),"-"))
 mores3<-grep(mores[3], colnames(TCGA1))
   if(length(mores3)!=0){
    mores2<-grep(mores[2], colnames(TCGA1)[mores3])
     if(length(mores2)!=0){
     moredata<-as.data.frame(TCGA1[,mores3]) 
     colnames(moredata)<-colnames(TCGA1)[mores3]
     moredataout<-cbind(moredataout,moredata)

}
}
}#39
dim(moredataout)

lessdataout<-TCGA1[,1]
for(i in 1:length(less)){
 lesss<-unlist(strsplit(as.character(less[i]),"-"))
 lesss3<-grep(lesss[3], colnames(TCGA1))
   if(length(lesss3)!=0){
    lesss2<-grep(lesss[2], colnames(TCGA1)[lesss3])
     if(length(lesss2)!=0){
     lessdata<-as.data.frame(TCGA1[,lesss3]) 
     colnames(lessdata)<-colnames(TCGA1)[lesss3]
     lessdataout<-cbind(lessdataout,lessdata)

}
}
}#49
dim(lessdataout)
GEO<-cbind(moredataout,lessdataout[,-1])

rownames(GEO)<-GEO[,1]
GEO<-GEO[,-1]

#GEO<-log2(GEO)######RNA-seq表达矩阵需要进行log2转化
###############创建设计矩阵##########
group1 <- rep('HIGH', (dim(moredataout)[2]-1))
group2 <- rep('LOW',(dim(lessdataout)[2]-1))
grouplist<-c(group1,group2)
table(grouplist)

design <- model.matrix(~0+factor(grouplist))#把分组变成因子形式
colnames(design)=levels(factor(grouplist))#改design列名为分组信息
exprSet=GEO#exprSet 为表达矩阵
rownames(design)=colnames(exprSet)#更改分组信息design行名为分组名，也就是样本名字

fit <- lmFit(exprSet, design)
cont.matrix <- makeContrasts("HIGH-LOW", levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
result<-na.omit(tempOutput)

result1<-cbind(rownames(result),result)
colnames(result1)<-c("protein",colnames(result))
GEO1<-cbind(rownames(GEO),GEO)
colnames(GEO1)<-c("gene",colnames(GEO))

final<-merge(GEO1, result1, by.x= "gene", by.y="protein")

write.table(final,"E:\\GBMdata\\GENE EXP\\RNA-seq\\result/HT_HG-U133A/HT_HG-U133A_gene_limma_result.txt",sep="\t",quote=FALSE, row.names = F)


library(gplots)
library(pheatmap)
library(RColorBrewer)
library(circlize)

#setwd("E:\\GBMdata\\GENE EXP\\RNA-seq\\result")

setwd("E:\\GBMdata\\GENE EXP\\RNA-seq\\result")
final<-read.table("HT_HG-U133A/HT_HG-U133A_gene_limma_result.txt",header=T,row.names=1,sep = "\t", quote = "")

up<-final[(final$P.Value<= 0.01 & final$logFC> 0.5),]
down<-final[(final$P.Value<= 0.01 & final$logFC< -0.5),]
gene2<-rbind(up,down)
for(i in 1:length(gene2)){
HIGHdata<-gene2[,1:(dim(moredataout)[2]-1)]
lowdata<-gene2[,(dim(moredataout)[2]):(dim(GEO)[2])]
#HIGHdata<-gene2[,1:23]
#lowdata<-gene2[,24:55]

dim(gene1)

gene1<-cbind(HIGHdata[,order(HIGHdata[i,])], lowdata[,order(lowdata[i,])])


test<-as.matrix(gene1[,1:(dim(GEO)[2])])
#colnames(test)<-rep(c("HIGH","LOW"),c((dim(moredataout)[2]-1),(dim(lessdataout)[2]-1)))
#test<-scale(test)
bk = unique(c(seq(-7, 0, length=100), seq(0, 7, length=100))) #scale范围
col = c(colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)) #颜色选取RdYlBu

#pheatmap还可以显示行或列的分组信息，支持多种分组；
annotation_col = data.frame(Group = factor(rep(c("HIGH","LOW"),c((dim(moredataout)[2]-1),(dim(lessdataout)[2]-1)))))
#annotation_col = data.frame(Group = factor(rep(c("HIGH","LOW"),c(23,32))))
rownames(annotation_col) =colnames(final)[1:88]
 
#还可以自己设定各个分组
ann_colors = list(Group = c(HIGH= "#1B9E77", LOW= "#D95F02"))
png(paste("E:\\GBMdata\\GENE EXP\\RNA-seq\\result\\HT_HG-U133A\\heatmap\\", i, ".png", sep=""))

pheatmap(test,
		annotation_col = annotation_col, #列分组
		cluster_rows = T,                #行聚类
		cluster_cols = F,                #列聚类
		annotation_colors = ann_colors,  #分组条颜色
		main="heatmap",                #标题
		#color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
		col =col,
		scale = "none"
		)
dev.off()

}


