mut<-read.table("E:\\GBMdata\\GENE EXP\\somatic mutation\\GBM_mc3_gene_level.txt",header=T,sep = "\t", quote = "", fill=T)

group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
HIGH<-group_name[group_name$group=="HIGH",1]
LOW<-group_name[group_name$group=="LOW",1]


HIGH_index_out<-c()
for(i in 1:length(HIGH)){
HIGH_index1<-grep(unlist(strsplit(as.character(HIGH[i]),"-"))[3],colnames(mut))
if(length(HIGH_index1)!=0){
HIGH_index_out<-c(HIGH_index_out, HIGH_index1)
}
}
length(HIGH_index_out)
HIGH_mut<-mut[,c(1,HIGH_index_out)]


LOW_index_out<-c()
for(i in 1:length(LOW)){
LOW_index1<-grep(unlist(strsplit(as.character(LOW[i]),"-"))[3],colnames(mut))
if(length(LOW_index1)!=0){
LOW_index_out<-c(LOW_index_out, LOW_index1)
}
}
length(LOW_index_out)
LOW_mut<-mut[,c(1,LOW_index_out)]
HIGH_mut_sum<-cbind(as.character(HIGH_mut[,1]), apply(HIGH_mut[,-1], 1,sum))
LOW_mut_sum<-cbind(as.character(LOW_mut[,1]), apply(LOW_mut[,-1], 1,sum))
colnames(HIGH_mut_sum)<-c("gene", "sum")
colnames(LOW_mut_sum)<-c("gene", "sum")

gene_position<-read.table("E:\\drug\\gene_position.txt",header=T,sep = "\t", quote = "")

mut_HIGH_position<-merge(gene_position,HIGH_mut_sum, by.x="DDX11L1", by.y="gene")[,c(2,3,1,6)]
mut_LOW_position<-merge(gene_position,LOW_mut_sum, by.x="DDX11L1", by.y="gene")[,c(2,3,1,6)]

