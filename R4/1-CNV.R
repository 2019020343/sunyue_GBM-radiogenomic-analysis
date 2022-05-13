rm(list=ls())
setwd("E:\\GBMdata\\GENE EXP\\CNV\\Gistic2_CopyNumber_Gistic2_all_data_by_genes")
cnv<-read.table("Gistic2_CopyNumber_Gistic2_all_data_by_gene.txt",header=T,sep = "\t", quote = "")
cnv_thresholded<-read.table("E:\\GBMdata\\GENE EXP\\CNV\\Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes\\Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",header=T,sep = "\t", quote = "")

group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
HIGH<-group_name[group_name$group=="HIGH",1]
LOW<-group_name[group_name$group=="LOW",1]


HIGH_index_out<-c()
for(i in 1:length(HIGH)){
HIGH_index1<-grep(unlist(strsplit(as.character(HIGH[i]),"-"))[3],colnames(cnv))
if(length(HIGH_index1)!=0){
HIGH_index_out<-c(HIGH_index_out, HIGH_index1)
}
}
length(HIGH_index_out)
HIGH_cnv<-cnv[,c(1,HIGH_index_out)]


LOW_index_out<-c()
for(i in 1:length(LOW)){
LOW_index1<-grep(unlist(strsplit(as.character(LOW[i]),"-"))[3],colnames(cnv))
if(length(LOW_index1)!=0){
LOW_index_out<-c(LOW_index_out, LOW_index1)
}
}
length(LOW_index_out)
LOW_cnv<-cnv[,c(1,LOW_index_out)]

cnv_out<-cbind(HIGH_cnv, LOW_cnv[,-1])

p_out<-c()
for(j in 1:nrow(cnv_out)){
	p<-wilcox.test(as.numeric(HIGH_cnv[j,-1]), as.numeric(LOW_cnv[j,-1]))[3]
	p_out<-c(p_out, p)
print(j)
}
fdr<-p.adjust(p_out,  method ="fdr")
cnv_out_fdr<-cbind(cnv_out, unlist(p_out), fdr )
colnames(cnv_out_fdr)<-c(colnames(cnv_out), "p", "fdr")
write.table(cnv_out_fdr,"E:\\GBMdata\\GENE EXP\\CNV\\cnv_out_fdr.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



