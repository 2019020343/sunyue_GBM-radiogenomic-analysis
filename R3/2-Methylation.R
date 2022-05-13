rm(list=ls())
library(gplots)
library(impute)
library(limma)

###  limma
setwd("E:\\GBMdata\\GENE EXP\\Methylation\\HumanMethylation450")
rt<-read.table("HumanMethylation450.txt",header=T,row.names=1,sep = "\t", quote = "", blank.lines.skip = F)
group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
less<-group_name[which(group_name$group=="LOW"),1]
more<-group_name[which(group_name$group=="HIGH"),1]
length(more)
length(less)

#LOW Group : sample match
moredataout<-rownames(rt)
for(i in 1:length(more)){
 mores<-unlist(strsplit(as.character(more[i]),"-"))
 mores3<-grep(mores[3], colnames(rt))
   if(length(mores3)!=0){
    mores2<-grep(mores[2], colnames(rt)[mores3])
     if(length(mores2)!=0){
     moredata<-as.data.frame(rt[,mores3]) 
     colnames(moredata)<-colnames(rt)[mores3]
     moredataout<-cbind(moredataout,moredata)

}
}
}#13
dim(moredataout)

#LOW Group : sample match
lessdataout<-rownames(rt)
for(i in 1:length(less)){
 lesss<-unlist(strsplit(as.character(less[i]),"-"))
 lesss3<-grep(lesss[3], colnames(rt))
   if(length(lesss3)!=0){
    lesss2<-grep(lesss[2], colnames(rt)[lesss3])
     if(length(lesss2)!=0){
     lessdata<-as.data.frame(rt[,lesss3]) 
     colnames(lessdata)<-colnames(rt)[lesss3]
     lessdataout<-cbind(lessdataout,lessdata)

}
}
}#8
dim(lessdataout)

## imoute  NA value 
methylation_data<-cbind(moredataout,lessdataout[,-1])
methylation_data<-as.matrix(methylation_data[,-1])
methylation_data<-impute.knn(methylation_data)
methylation_data<-as.data.frame(methylation_data$data)

## DMA 
out<-c()
for(i in 1:nrow(methylation_data)){
	moremethy<-methylation_data[i,1:(ncol(moredataout)-1)]
	lessmethy<-methylation_data[i,ncol(moredataout):ncol(methylation_data)]
      p<-t.test(as.numeric(moremethy),as.numeric(lessmethy))[[3]]
      diff_moreCless<-mean(as.numeric(moremethy))-mean(as.numeric(lessmethy))
      outdata<-cbind(p, diff_moreCless)
	out<-rbind(outdata,out)
print(i)

}
final<-cbind(methylation_data, out)
final1<-cbind(rownames(final), final)
colnames(final1)<-c("ID", colnames(final))
multi<-read.table("multi.txt",header=T,sep = "\t", quote = "", blank.lines.skip = F)#multi==0
name1 <- multi[which(multi[,2]!=0),1]

snpAchr<-read.table("snpAchr.txt",header=T,sep = "\t", quote = "", check.names=F)##snpAchr : SNP_COMMON=F chrX Y M=F 
name2 <- snpAchr[which(snpAchr[,3]=="TRUE" | snpAchr[,2]=="chrX" | snpAchr[,2]=="chrY" | snpAchr[,2]=="chrM"),1]
name3 <- unique(c(as.character(name1),as.character(name2)))

A<-setdiff(final1[,1],  name3)

data2 <- final1[A,]
dara3<-p.adjust(data2$p, method = "fdr")
dara5<-data.frame(data2, dara3)
dara4<-dara5[which(dara5$dara3<0.05 & abs(dara5$diff_moreCless)>0.2), ]
write.table(dara5,"450k_fdr_all.txt",sep="\t",quote=FALSE, row.name=F, col.name=T)


#####  gene  match 
setwd("E:\\GBMdata\\GENE EXP\\Methylation\\HumanMethylation450")
cgchr<-read.table("illuminaMethyl450_hg19_GPL16304_TCGAlegacy.txt",header=F,sep = "\t", quote = "")
merhy<-read.table("450k_fdr_all.txt",header=T,sep = "\t", quote = "")
colnames(cgchr)<-c("cgID", "genesymbol", "chr", "start", "end", "state")
upcg<-merhy[(merhy$p<= 0.05 & merhy$diff_moreCless>= 0.2),]
downcg<-merhy[(merhy$p<= 0.05 & merhy$diff_moreCless<= -0.2),]
me2<-rbind(upcg,downcg)

gene_me<-merge(cgchr,me2,by.x="cgID", by.y="ID")
head(gene_me)
gene_me_out<-c()
for(i in 1:nrow(gene_me)){
	genesymbol1<-unlist(strsplit(as.character(gene_me$genesymbol[i]),","))
	if(length(genesymbol1)>=1){
		gene_me_1<-cbind(genesymbol1, gene_me[i,])
}
gene_me_out<-rbind(gene_me_out,gene_me_1)
}
write.table(gene_me_out,"cg_gene_match_result.txt",sep="\t",quote=FALSE, row.name=F, col.name=T)





################ TI cg
rm(list=ls())
setwd("E:\\GBMdata\\GENE EXP\\Methylation\\HumanMethylation450")
LOW_LOW<-read.table("cg_gene_match_result.txt",header=T,sep = "\t", quote = "")
hkGENE<-read.table("E:\\immune\\2MGN\\hkGENE18.txt",header=T,sep = "\t", quote = "")
TIgene_LOW<-read.table("E:\\drug\\TIGENE_LOWname_21.txt",header=T,sep = "\t", quote = "")

LOW_hyper<-LOW_LOW[HIGH_LOW$p<0.05 & HIGH_LOW$diff_moreCless>0.2, ]
HIGH_hyper[,1]
genenames<-unique(HIGH_hyper[,1])  ##取向量test中的重复值
genenames<-intersect(hkGENE[,1],HIGH_hyper[,1])
sum(is.na(match(hkGENE[,1], HIGH_hyper[,1]))==0)

colnames(HIGH_hyper)
HIGH_hyper_1<-c()
for(k in 1:length(genenames)){
	cg1<-HIGH_hyper[which(HIGH_hyper[,1]==genenames[k]),c(1:2,8:20)]
	HIGH_hyper_1<-rbind(HIGH_hyper_1, cg1)

}
HIGH_hyper<-unique(HIGH_hyper_1)
head(HIGH_hyper)

LOW_hyper<-HIGH_LOW[HIGH_LOW$p<0.05 & HIGH_LOW$diff_moreCless<(-0.2), ]
genenames_LOW<-unique(LOW_hyper[,1])
genenames_LOW<-intersect(hkGENE[,1],LOW_hyper[,1])
LOW_hyper_1<-c()
for(k in 1:length(genenames_LOW)){
	LOW_cg1<-LOW_hyper[which(LOW_hyper[,1]==genenames_LOW[k]), c(1:2,21:28)]
	LOW_hyper_1<-rbind(LOW_hyper_1, LOW_cg1)

}
LOW_hyper<-unique(LOW_hyper_1)


write.table(HIGH_hyper,"E:\\drug\\GDSC\\methylation\\hyper_HIGH.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(LOW_hyper,"E:\\drug\\GDSC\\methylation\\hyper_LOW.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



## CELL SAMPLE:MATCH
rm(list=ls())
HIGH_hyper<-read.table("E:\\drug\\GDSC\\methylation\\hyper_HIGH.txt",header=T,sep = "\t", quote = "")
LOW_hyper<-read.table("E:\\drug\\GDSC\\methylation\\hyper_LOW.txt",header=T,sep = "\t", quote = "")

GBM_CELLnames<-read.table("E:\\drug\\GDSC\\cellline_names.txt",header=T,sep = "\t", quote = "")
colname<-read.table("D:\\additional data\\methy\\colname.txt",header=F,sep = "\t", quote = "")
colname<-colname[1,1:1029] #细胞系的beta值
colname_1<-as.data.frame(unlist(sapply(colname, function(x){unlist(strsplit(as.character(x),"_"))[1]})))#细胞系的名字
myfun<-function(x,y){which(y[,1]==as.character(x))}
a<-as.data.frame(cbind(as.character(GBM_CELLnames[,1]),sapply(GBM_CELLnames[,1], myfun, colname_1)))
a<-a[-5,]
HIGH_cg_all_out<-c()
LOW_cg_all_out<-c()
for(h in 4:17){
	GDSC_cg_data<-read.table(paste("D:\\additional data\\methy\\",h-1,".txt",sep=""),header=F,sep = "\t", quote = "")
	GDSC_cg_GBM<-GDSC_cg_data[,c(2,unlist(a[,2]))]
	HIGH_cg_all_1<-merge(HIGH_hyper, GDSC_cg_GBM, by.x="cgID", by.y="V2")	
	LOW_cg_all_1<-merge(LOW_hyper, GDSC_cg_GBM, by.x="cgID", by.y="V2")	
	HIGH_cg_all_out<-rbind(HIGH_cg_all_out, HIGH_cg_all_1)
	LOW_cg_all_out<-rbind(LOW_cg_all_out, LOW_cg_all_1)


print(h)
}
colnames(HIGH_cg_all_out)<-c(colnames(HIGH_cg_all_out)[1:15],unlist(a[,1]))
colnames(LOW_cg_all_out)<-c(colnames(LOW_cg_all_out)[1:10],unlist(a[,1]))

head(TIGENE_LOW_celllines_GBM)
#筛选GBM细胞系 和 样本 cg的beta 
TIGENE_HIGH_samples_GBM<-HIGH_cg_all_out[,1:15]
TIGENE_HIGH_celllines_GBM<-HIGH_cg_all_out[,c(1:2,16:ncol(HIGH_cg_all_out))]

TIGENE_LOW_samples_GBM<-LOW_cg_all_out[,1:10]
TIGENE_LOW_celllines_GBM<-LOW_cg_all_out[,c(1:2,11:ncol(LOW_cg_all_out))]

TIGENE_HIGH_samples_GBM_scale<-cbind(TIGENE_HIGH_samples_GBM[,1:2], scale(TIGENE_HIGH_samples_GBM[,-(1:2)],center=T,scale=T))
TIGENE_HIGH_celllines_GBM_scale<-cbind(TIGENE_HIGH_celllines_GBM[,1:2], scale(TIGENE_HIGH_celllines_GBM[,-(1:2)],center=T,scale=T))


TIGENE_LOW_samples_GBM_scale<-cbind(TIGENE_LOW_samples_GBM[,1:2], scale(TIGENE_LOW_samples_GBM[,-(1:2)],center=T,scale=T))
TIGENE_LOW_celllines_GBM_scale<-cbind(TIGENE_LOW_celllines_GBM[,1:2], scale(TIGENE_LOW_celllines_GBM[,-(1:2)],center=T,scale=T))


#皮尔森相关系数r，相异性矩阵1-r
samples_celllines_HIGH<-c()
max<-matrix(100,ncol=ncol(TIGENE_HIGH_samples_GBM), nrow=ncol(TIGENE_HIGH_celllines_GBM))
max[1,]<-colnames(TIGENE_HIGH_samples_GBM)
max[,1]<-colnames(TIGENE_HIGH_celllines_GBM)
for(i in 3:ncol(TIGENE_HIGH_samples_GBM)){
	for(j in 3:ncol(TIGENE_HIGH_celllines_GBM)){
	corRESULT<-cor.test(as.numeric(TIGENE_HIGH_samples_GBM_scale[,i]),as.numeric(TIGENE_HIGH_celllines_GBM_scale[,j]),method="pearson")
	p_value<-corRESULT[3]
	r_value<-corRESULT[4]
	max[j,i]<-1-as.numeric(r_value)
	if(r_value>0.8 & p_value<0.01){
samples_celllines_HIGH_1<-cbind(colnames(TIGENE_HIGH_samples_GBM)[i],colnames(TIGENE_HIGH_celllines_GBM)[j],r_value)
samples_celllines_HIGH<-rbind(samples_celllines_HIGH,samples_celllines_HIGH_1)
}	
}	
}

#皮尔森相关系数r，相异性矩阵1-r
samples_celllines_LOW<-c()
max_LOW<-matrix(100,ncol=ncol(TIGENE_LOW_samples_GBM), nrow=ncol(TIGENE_LOW_celllines_GBM))
max_LOW[1,]<-colnames(TIGENE_LOW_samples_GBM)
max_LOW[,1]<-colnames(TIGENE_LOW_celllines_GBM)
for(i in 3:ncol(TIGENE_LOW_samples_GBM)){
	for(j in 3:ncol(TIGENE_LOW_celllines_GBM)){
	corRESULT<-cor.test(as.numeric(TIGENE_LOW_samples_GBM_scale[,i]),as.numeric(TIGENE_LOW_celllines_GBM_scale[,j]),method="pearson")
	p_value<-corRESULT[3]
	r_value<-corRESULT[4]
	max_LOW[j,i]<-1-as.numeric(r_value)
	if(r_value>0.8 & p_value<0.01){
samples_celllines_LOW_1<-cbind(colnames(TIGENE_LOW_samples_GBM)[i],colnames(TIGENE_LOW_celllines_GBM)[j],r_value)
samples_celllines_LOW<-rbind(samples_celllines_LOW,samples_celllines_LOW_1)
}	
}	
}
samples_celllines_HIGH_methy<-c()
for(i in 1:nrow(samples_celllines_HIGH)){
	a1<-unlist(samples_celllines_HIGH[i,])
	samples_celllines_HIGH_methy<-rbind(samples_celllines_HIGH_methy, a1)
}
samples_celllines_LOW_methy<-c()
for(i in 1:nrow(samples_celllines_LOW)){
	a1<-unlist(samples_celllines_LOW[i,])
	samples_celllines_LOW_methy<-rbind(samples_celllines_LOW_methy, a1)
}
write.table(samples_celllines_HIGH_methy,"E:\\drug\\samples_celllines_HIGH_methy.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
write.table(samples_celllines_LOW_methy,"E:\\drug\\samples_celllines_LOW_methy.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
intersect(samples_celllines_HIGH_methy[,2], samples_celllines_HIGH[,2])
