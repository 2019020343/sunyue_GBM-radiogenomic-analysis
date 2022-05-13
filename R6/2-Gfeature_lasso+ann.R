################# gene ###############
rm(list=ls())
library(e1071)
library(pROC)
library(plyr)
library(glmnet)

group<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",sep="\t",header=T)
LGG_group<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA\\group.txt",sep="\t",header=T, row.names = NULL)
### 21  genes
TIgene_HIGH<-read.table("E:\\drug\\TIGENE_HIGHname_21.txt",header=T,sep = "\t", quote = "", check.names=F)
HIGH_LOW<-read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\result/HT_HG-U133A/HT_HG-U133A_gene_limma_result.txt",header=T,sep = "\t", quote = "")
LGG_HIGH_LOW<-read.table("E:\\LGGdata\\HiSeqV2.txt",header=T,sep = "\t", quote = "")


TIgene_HIGH_exp<-merge(TIgene_HIGH, HIGH_LOW, by.x="symbol", by.y="gene")[,c(1,6:93)]
TIgene_HIGH_exp<-t(TIgene_HIGH_exp)

LGG_HIGH_LOW<-LGG_HIGH_LOW[,c(1,grep("01$", colnames(LGG_HIGH_LOW)))]
LGG_TIgene_HIGH_exp<-merge(TIgene_HIGH, LGG_HIGH_LOW, by.x="symbol", by.y="sample")[,c(1,6:521)]
LGG_TIgene_HIGH_exp<-t(LGG_TIgene_HIGH_exp)



myfun <- function(x){
  unlist(strsplit(x,".",fixed = T))[3]
}
samplename<-sapply(rownames(TIgene_HIGH_exp),myfun)
LGG_samplename<-sapply(rownames(LGG_TIgene_HIGH_exp),myfun)
TIgene_HIGH_exp_1<-as.data.frame(TIgene_HIGH_exp)
LGG_TIgene_HIGH_exp_1<-as.data.frame(LGG_TIgene_HIGH_exp)

#write.table(TIgene_HIGH_exp_1,"E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\TIgene_HIGH_exp_1.txt",sep="\t",quote=FALSE, row.name=T, col.name=F)
#TIgene_HIGH_exp_1<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\TIgene_HIGH_exp_1.txt",sep="\t",header=T)


for (h in 2:length(samplename)) {
  label<-group$group[grep(samplename[h], group[,1])]
  if(length(label)!=0){
    TIgene_HIGH_exp_1[h, 22] = as.character(label)
  }
}
TIgene_HIGH_exp_1<-TIgene_HIGH_exp_1[-1,]

LGG_TIgene_HIGH_exp_out<-c()
for (h in 2:length(LGG_samplename)) {
  label<-LGG_group$group[grep(LGG_samplename[h], LGG_group[,1])]
  if(length(label)!=0){
    LGG_TIgene_HIGH_exp_2<-cbind(LGG_TIgene_HIGH_exp_1[h, ], as.character(label))
    LGG_TIgene_HIGH_exp_out<-rbind(LGG_TIgene_HIGH_exp_out, LGG_TIgene_HIGH_exp_2)
  }
}
LGG_TIgene_HIGH_exp_1<-LGG_TIgene_HIGH_exp_out



##### ±ê×¼»¯  #######
data_MF_exp_scale_1<-as.numeric(as.matrix(TIgene_HIGH_exp_1[,-ncol(TIgene_HIGH_exp_1)]))
data_MF_exp_scale_2<-matrix(data_MF_exp_scale_1, ncol = ncol(TIgene_HIGH_exp_1)-1)
data_MF_exp_scale<-scale(data_MF_exp_scale_2 , center=T,scale=T)

data_MF_exp<-as.data.frame(cbind(data_MF_exp_scale, as.character(TIgene_HIGH_exp_1[,ncol(TIgene_HIGH_exp_1)]))) 
colnames(data_MF_exp)<-c(as.character(TIgene_HIGH$symbol), "group")
rownames(data_MF_exp)<-rownames(TIgene_HIGH_exp_1)
write.table(data_MF_exp,"E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\data_MF_exp.txt",sep="\t",quote=FALSE, row.name=T, col.name=T)


LGG_data_MF_exp_scale_1<-as.numeric(as.matrix(LGG_TIgene_HIGH_exp_1[,-ncol(LGG_TIgene_HIGH_exp_1)]))
LGG_data_MF_exp_scale_2<-matrix(LGG_data_MF_exp_scale_1, ncol = ncol(LGG_TIgene_HIGH_exp_1)-1)
LGG_data_MF_exp_scale<-scale(LGG_data_MF_exp_scale_2 , center=T,scale=T)

LGG_data_MF_exp<-as.data.frame(cbind(LGG_data_MF_exp_scale, as.character(LGG_TIgene_HIGH_exp_1[,ncol(LGG_TIgene_HIGH_exp_1)]))) 
colnames(LGG_data_MF_exp)<-c(as.character(TIgene_HIGH$symbol), "group")
rownames(LGG_data_MF_exp)<-rownames(LGG_TIgene_HIGH_exp_1)
write.table(LGG_data_MF_exp,"E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\lasso\\LGG_data_MF_exp.txt",sep="\t",quote=FALSE, row.name=T, col.name=T)






##### model  #########

out_08<-c()
out_8<-c()
seedset<-sample(1000:2500000, 1000, replace = F)
for (j in 1:1000) {
  set.seed(seedset[j])
  n<-length(TIgene_HIGH_exp_1[,1])
  samp<-sample(1:n,n/5)
  
  GBM_train<-data_MF_exp[-samp, ]
  GBM_test<-data_MF_exp[samp, ]

  
  x1<-GBM_test[,1:21]
  x1[sapply(x1, is.factor)] <- lapply(x1[sapply(x1, is.factor)], function(m) as.numeric(as.character(m)))
  y=GBM_test[,22]
  GBM_test<-cbind(x1,y)
  
  
  x=GBM_train[,1:21]
  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))
  y=GBM_train[,22]
  GBM_train<-cbind(x,y)
  
  x=LGG_data_MF_exp[,1:21]
  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))
  y=LGG_data_MF_exp[,22]
  LGG_anntest<-cbind(x,y)
  
  
    ir.nn2<-nnet(y~.,data=GBM_train,  linout=F,size=5,decay=5e-4,maxit=100,trace=F,rang=0.7)
    
    pre_ann<-predict(ir.nn2, GBM_test, type = "class")
    LGG_pre_ann<-predict(ir.nn2, LGG_anntest, type = "class")
    
    
    pre_ann<-ifelse(pre_ann=="HIGH", 1, 2)
    GBM_test$y<-ifelse(GBM_test$y=="HIGH", 1, 2)
    
    LGG_pre_ann<-ifelse(LGG_pre_ann=="HIGH", 1, 2)
    LGG_anntest$y<-ifelse(LGG_anntest$y=="HIGH", 1, 2)
    
    GBM_auc <- roc(GBM_test$y,pre_ann )$auc
    LGG_auc <- roc(LGG_anntest$y,LGG_pre_ann)$auc
    
    out_1<-cbind(j,seedset[j],GBM_auc, LGG_auc)
    out_08<-rbind(out_08, out_1)
    print(j)
    
}

###########  final  model ########
#rm(list=ls())
data_MF_exp<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\data_MF_exp.txt",header=T,sep = "\t", quote = "", check.names=F)
LGG_data_MF_exp<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\lasso\\LGG_data_MF_exp.txt",header=T,sep = "\t", quote = "", check.names=F)



set.seed(1223965)
n<-length(data_MF_exp[,1])
samp<-sample(1:n,n/5)
GBM_train<-data_MF_exp[-samp, ]
GBM_test<-data_MF_exp[samp, ]

x1<-GBM_test[,1:21]
x1[sapply(x1, is.factor)] <- lapply(x1[sapply(x1, is.factor)], function(m) as.numeric(as.character(m)))
y=GBM_test[,22]
GBM_test_G<-cbind(x1,y)


x=GBM_train[,1:21]
x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))
y=GBM_train[,22]
GBM_train<-cbind(x,y)

x=LGG_data_MF_exp[,1:21]
x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))
y=LGG_data_MF_exp[,22]
LGG_anntest_G<-cbind(x,y)


ir.nn2<-nnet(y~.,data=GBM_train,  linout=F,size=5,decay=5e-4,maxit=100,trace=F, rang=0.7)

pre_ann_G<-predict(ir.nn2, GBM_test_G, type = "class")
LGG_pre_ann_G<-predict(ir.nn2, LGG_anntest_G, type = "class")


pre_ann_G<-ifelse(pre_ann_G=="HIGH", 1, 2)
GBM_test_G$y<-ifelse(GBM_test_G$y=="HIGH", 1, 2)

LGG_pre_ann_G<-ifelse(LGG_pre_ann_G=="HIGH", 1, 2)
LGG_anntest_G$y<-ifelse(LGG_anntest_G$y=="HIGH", 1, 2)

GBM_auc_G <- roc(GBM_test_G$y,pre_ann_G )
LGG_auc_G <- roc(LGG_anntest_G$y,LGG_pre_ann_G)







  
  

