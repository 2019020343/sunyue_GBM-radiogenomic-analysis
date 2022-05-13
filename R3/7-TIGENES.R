rm(list=ls())
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(RColorBrewer)
##############   gene  exp
HIGH_LOW<-read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\result/HT_HG-U133A/HT_HG-U133A_gene_limma_result.txt",header=T,sep = "\t", quote = "")
hkGENE<-read.table("E:\\immune\\2MGN\\hkGENE18.txt",header=T,sep = "\t", quote = "")
HIGH_down<-HIGH_LOW[HIGH_LOW$P.Value<0.01 & HIGH_LOW$logFC<0, ]
LOW_down<-HIGH_LOW[HIGH_LOW$P.Value<0.01 & HIGH_LOW$logFC>0, ]
TIgene_HIGH<-intersect(HIGH_down[,1],hkGENE[,1])
TIgene_LOW<-intersect(LOW_down[,1],hkGENE[,1])
dim(HIGH_down)
dim(LOW_down)

length(TIgene_HIGH)
length(TIgene_LOW)
#write.table(TIgene_HIGH,"E:\\drug\\TIGENE_HIGHname_21.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

#手动筛选中21个基因的ENSG ID
TIgene_HIGH<-read.table("E:\\drug\\TIGENE_HIGHname_21.txt",header=T,sep = "\t", quote = "", check.names=F)

TIGENE_HIGH_celllines<-read.table("E:\\drug\\GDSC\\sanger1018_brainarray_ensemblgene_rma.txt",header=T,sep = "\t", quote = "", check.names=F, fill=T)
cellline_names<-read.table("E:\\drug\\GDSC\\cellline_names.txt",header=T,sep = "\t", quote = "", check.names=F, fill=T)


TIGENE_HIGH_celllines_1<-merge(TIGENE_HIGH_celllines, TIgene_HIGH, by.x="ensembl_gene", by.y="ENSG")
dim(TIGENE_HIGH_celllines_1)



#TIGENE_LOW_celllines<-read.table("E:\\drug\\TIGENE_LOW_4.txt",header=T,sep = "\t", quote = "")
#筛选中枢神经系统细胞系即GBM或LGG细胞系 TIGENE的exp
TIGENE_HIGH_celllines_GBM<-TIGENE_HIGH_celllines_1[,c(ncol(TIGENE_HIGH_celllines_1),unlist(sapply(cellline_names[,2],myfun<-function(x){which(colnames(TIGENE_HIGH_celllines_1)==x)})))]
TIGENE_HIGH_celllines_GBM<-TIGENE_HIGH_celllines_GBM[order(TIGENE_HIGH_celllines_GBM[,1]),]
#TIGENE_LOW_celllines_GBM<-TIGENE_LOW_celllines[,c(2,grep("CENTRAL_NERVOUS_SYSTEM",colnames(TIGENE_LOW_celllines)))]

dim(TIGENE_HIGH_celllines_GBM)
#筛选GBM样本 TIGENE的exp
TIGENE_HIGH_samples_GBM<-merge(HIGH_LOW, TIgene_HIGH, by.x="gene",by.y="symbol")[,1:40]
TIGENE_HIGH_samples_GBM<-TIGENE_HIGH_samples_GBM[order(TIGENE_HIGH_samples_GBM[,1]), ]
#TIGENE_LOW_samples_GBM<-merge(as.data.frame(TIgene_LOW), HIGH_LOW,by.x="TIgene_LOW",by.y="gene")[,c(1,25:56)]


#标准化
TIGENE_HIGH_samples_GBM_scale<-cbind(TIGENE_HIGH_samples_GBM[,1], scale(TIGENE_HIGH_samples_GBM[,-1],center=T,scale=T))
TIGENE_HIGH_celllines_GBM_scale<-cbind(TIGENE_HIGH_celllines_GBM[,1], scale(TIGENE_HIGH_celllines_GBM[,-1],center=T,scale=T))
dim(TIGENE_HIGH_celllines_GBM)
#皮尔森相关系数r，相异性矩阵1-r
samples_celllines_HIGH<-c()
samples_celllines_unHIGH<-c()

max<-matrix(100,ncol=ncol(TIGENE_HIGH_samples_GBM_scale), nrow=ncol(TIGENE_HIGH_celllines_GBM_scale))
max[1,]<-colnames(TIGENE_HIGH_samples_GBM_scale)
max[,1]<-colnames(TIGENE_HIGH_celllines_GBM_scale)
for(i in 2:ncol(TIGENE_HIGH_samples_GBM_scale)){
	for(j in 2:ncol(TIGENE_HIGH_celllines_GBM_scale)){
	r_value<-cor(TIGENE_HIGH_samples_GBM_scale[,i],TIGENE_HIGH_celllines_GBM_scale[,j])
	max[j,i]<-as.numeric(r_value)
	if(r_value>0.8){
clname<-as.character(cellline_names[which(cellline_names[,2]==colnames(TIGENE_HIGH_celllines_GBM_scale)[j]),1])
samples_celllines_HIGH_1<-cbind(colnames(TIGENE_HIGH_samples_GBM_scale)[i],clname,r_value)
samples_celllines_HIGH<-rbind(samples_celllines_HIGH,samples_celllines_HIGH_1)
}else{
clname<-as.character(cellline_names[which(cellline_names[,2]==colnames(TIGENE_HIGH_celllines_GBM_scale)[j]),1])
samples_celllines_unHIGH_1<-cbind(colnames(TIGENE_HIGH_samples_GBM_scale)[i],clname,r_value)
samples_celllines_unHIGH<-rbind(samples_celllines_unHIGH,samples_celllines_unHIGH_1)
}	
}	
}
#相关系>0.8的样本-细胞系对
r8_10<-unlist(unique(samples_celllines_HIGH[,2]))
r0_8<-setdiff(cellline_names[,1], r8_10 )
colnames(samples_celllines_HIGH)<-c("samples", "celllines", "r_value")

#write.table(samples_celllines_HIGH, "E:\\drug\\samples_celllines_HIGH.txt", row.name=F, sep="\t",quote=FALSE)

##################  heatmap  ###############################
r8_10<-as.data.frame(r8_10)
colnames(r8_10)<-"Sample Name"
r8_10_name<-merge(r8_10, cellline_names, by.x="Sample Name", by.y="Sample Name")
indexout<-c()
for(i in 1:nrow(r8_10_name)){
	index<-which(colnames(TIGENE_HIGH_celllines_GBM)==r8_10_name[i,2])
	indexout<-c(indexout, index)
}

samples= as.matrix(TIGENE_HIGH_samples_GBM[,-1])
samples<-scale(samples)
rownames(samples)<-TIGENE_HIGH_samples_GBM[,1]


celllines= as.matrix(TIGENE_HIGH_celllines_GBM[,-1])
celllines<-scale(celllines)
rownames(celllines)<-TIGENE_HIGH_celllines_GBM[,1]


celllines_1= as.matrix(TIGENE_HIGH_celllines_GBM[,indexout])
celllines_1<-scale(celllines_1)
rownames(celllines_1)<-TIGENE_HIGH_celllines_GBM[,1]


Heatmap(samples, name = "samples", col=colorRampPalette( c("#6E98D6", "white", "#FBBD8E"))(200),
show_row_names = FALSE, show_column_names = F) +

#BEA0F9  ABBE8C
Heatmap(celllines_1, name = "celllines", col=colorRampPalette(c("#6E98D6", "white", "#FBBD8E"))(200), 
        column_names_gp = gpar(fontsize = 8),column_names_rot = 90,show_column_names = F)

##################### circlize ##############
library(dplyr)
library(statnet)
library(circlize)

data<-as.data.frame(rbind(samples_celllines_HIGH, samples_celllines_unHIGH))
data1<-data[which(as.numeric(as.matrix(data[,3]))> 0.8),]

chordDiagram(data1, 
             annotationTrack = c('grid','axis'),
             directional = TRUE,
             diffHeight = 0,
             #grid.col = grid.col, 
             transparency = 0.5,                       
             preAllocateTracks = list(
               track.height = uh(4, "mm"),
               track.margin = c(uh(4, "mm"), 0))
             )


