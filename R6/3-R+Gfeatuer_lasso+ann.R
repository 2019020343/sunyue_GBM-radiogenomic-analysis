rm(list=ls())
GBM_data_MF<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\0.73\\F_M_matrix_BZ.txt",sep="\t",header=T, row.names = 1)
LGG_data_MF<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\lasso\\F_M_matrix_BZ.txt",sep="\t",header=T, row.names = 1)

### 21  genes
TIgene_HIGH<-read.table("E:\\drug\\TIGENE_HIGHname_21.txt",header=T,sep = "\t", quote = "", check.names=F)

HIGH_LOW<-read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\result/HT_HG-U133A/HT_HG-U133A_gene_limma_result.txt",header=T,sep = "\t", quote = "")
HIGH_LOW2<-scale(HIGH_LOW[,-1])     
HIGH_LOW3<-cbind(as.character(HIGH_LOW[,1]), HIGH_LOW2)

LGG_HIGH_LOW<-read.table("E:\\LGGdata\\HiSeqV2.txt",header=T,sep = "\t", quote = "")
LGG_HIGH_LOW<-LGG_HIGH_LOW[,c(1,grep("01$", colnames(LGG_HIGH_LOW)))]
LGG_HIGH_LOW2<-scale(LGG_HIGH_LOW[,-1])     
LGG_HIGH_LOW3<-cbind(as.character(LGG_HIGH_LOW[,1]), LGG_HIGH_LOW2)

TIgene_HIGH_exp<-merge(TIgene_HIGH, HIGH_LOW3, by.x="symbol", by.y="V1")[,c(1,6:93)]
LGG_TIgene_HIGH_exp<-merge(TIgene_HIGH, LGG_HIGH_LOW3, by.x="symbol", by.y="V1")[,c(1,6:521)]

data_MF_exp_1<-c()
for (k in 1:nrow(GBM_data_MF)) {
  samples1<-strsplit(as.character(rownames(GBM_data_MF)[k]), "_")[[1]][1]
  samples2<-strsplit(as.character(samples1), "-")[[1]][2]
  geneexp_match<-TIgene_HIGH_exp[,grep(samples2, colnames(TIgene_HIGH_exp))]
  if(length(geneexp_match)!=0){
    data_MF_exp1<-unlist(t(c(rownames(GBM_data_MF)[k], GBM_data_MF[k,-ncol(GBM_data_MF)], as.matrix(geneexp_match), as.character(GBM_data_MF[k,ncol(GBM_data_MF)]))))
    data_MF_exp_1<-rbind(data_MF_exp_1, data_MF_exp1)
  }  
}

LGG_data_MF_exp_1<-c()
for (k in 1:nrow(LGG_data_MF)) {
  samples1<-strsplit(as.character(rownames(LGG_data_MF)[k]), "_")[[1]][1]
  samples2<-strsplit(as.character(samples1), "-")[[1]][2]
  geneexp_match<-LGG_TIgene_HIGH_exp[,grep(samples2, colnames(LGG_TIgene_HIGH_exp))]
  if(length(geneexp_match)!=0){
    LGG_data_MF_exp1<-unlist(t(c(rownames(LGG_data_MF)[k], LGG_data_MF[k,-ncol(LGG_data_MF)],as.matrix(geneexp_match),as.character(LGG_data_MF[k,ncol(LGG_data_MF)]))))
    LGG_data_MF_exp_1<-rbind(LGG_data_MF_exp_1, LGG_data_MF_exp1)
  }  
}

rownames(data_MF_exp_1)<-data_MF_exp_1[,1]
data_MF_exp_1<-as.data.frame(data_MF_exp_1[,-1]) 
colnames(data_MF_exp_1)<-c(colnames(GBM_data_MF)[-ncol(GBM_data_MF)], as.character(TIgene_HIGH_exp[,1]), "group")
  
rownames(LGG_data_MF_exp_1)<-LGG_data_MF_exp_1[,1]
LGG_data_MF_exp_1<-as.data.frame(LGG_data_MF_exp_1[,-1])
colnames(LGG_data_MF_exp_1)<-c(colnames(LGG_data_MF)[-ncol(LGG_data_MF)], as.character(LGG_TIgene_HIGH_exp[,1]), "group")

#write.table(data_MF_exp_1,"E:\\Mask_RCNN_master\\pyradiomics\\test\\Lasso\\RG_matrix_BZ.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)
#write.table(LGG_data_MF_exp_1,"E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\Lasso\\RG_matrix_BZ.txt",col.names = T, row.names = T,sep = "\t" ,append = FALSE, quote = F)
###################################   ann   #################
#rm(list=ls())
library(nnet)
library(glmnet)
library(pROC)
library(plyr)
library(glmnet)
data_MF_exp_1<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\Lasso\\RG_matrix_BZ.txt",header=T,sep = "\t", quote = "", check.names=F)
LGG_data_MF_exp_1<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\Lasso\\RG_matrix_BZ.txt",header=T,sep = "\t", quote = "", check.names=F)

out_08<-c()
out_8<-c()
seedset<-sample(1000:2500000, 500, replace = F)
for (j in 1:500) {
  set.seed(seedset[j])
  n<-length(data_MF_exp_1[,1])
  samp<-sample(1:n,n/5)
  GBM_train<-data_MF_exp_1[-samp, ]
  GBM_test<-data_MF_exp_1[samp, ]

  x=GBM_train[,1:107]
  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))
  
  y=as.numeric(GBM_train[,108])
  cv.fit<-cv.glmnet(data.matrix(x) ,y,family="binomial") 
  plot(cv.fit)
  
  fit<-glmnet(data.matrix(x),y,family="binomial", nlambda=500, alpha=1)
  # print(fit)
  # plot(fit)
  plot(fit, xvar="lambda", label=TRUE)
  coefficients<-unlist(coef(fit,s=cv.fit$lambda.min))
  Active.Index<-which(coefficients!=0) #系数不为0的特征索引
  Active.feature<-rownames(coefficients)[which(coefficients!=0)]
  Active.coefficients<-coefficients[Active.Index] #系数不为0的特征系数值
  if(length(Active.feature)>8){
    anntrain<-GBM_train[,c(Active.Index-1, ncol(data_MF_exp_1))]
    anntest<-GBM_test[,c(Active.Index-1, ncol(data_MF_exp_1))]
    LGG_anntest<-LGG_data_MF_exp_1[,c(Active.Index-1, ncol(LGG_data_MF_exp_1))]
    
    
    ir.nn2<-nnet(group~., data = anntrain,linout=F,size=10,decay=0.001,maxit=100,trace=F, MaxNWts = 100000)
    pre_ann<-predict(ir.nn2, anntest, type = "class")
    LGG_pre_ann<-predict(ir.nn2, LGG_anntest, type = "class")
    
    
    pre_ann<-ifelse(pre_ann=="HIGH", 1, 0)
    anntest$group<-ifelse(anntest$group=="HIGH", 1, 0)
    
    LGG_pre_ann<-ifelse(LGG_pre_ann=="HIGH", 1, 0)
    LGG_anntest$group<-ifelse(LGG_anntest$group=="HIGH", 1, 0)
    
    GBM_auc <- roc(anntest$group,pre_ann )$auc
    LGG_auc <- roc(LGG_anntest$group,LGG_pre_ann)$auc
    
    if(GBM_auc>0.8 & LGG_auc>0.7){
      out_1<-cbind(j,seedset[j],GBM_auc, LGG_auc)
      print(j)
    }
    out_1<-cbind(j,seedset[j],GBM_auc, LGG_auc)
    out_08<-rbind(out_08, out_1)

    
  }
  
}

###################  final  MODEL#############
#rm(list=ls())
library(nnet)
library(glmnet)
library(pROC)
library(plyr)
library(glmnet)
data_MF_exp_1<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\Lasso\\RG_matrix_BZ.txt",header=T,sep = "\t", quote = "", check.names=F)
LGG_data_MF_exp_1<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\Lasso\\RG_matrix_BZ.txt",header=T,sep = "\t", quote = "", check.names=F)

set.seed(1228917)
n<-length(data_MF_exp_1[,1])
samp<-sample(1:n,n/5)
GBM_train<-data_MF_exp_1[-samp, ]
GBM_test<-data_MF_exp_1[samp, ]

x=GBM_train[,1:107]
x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))

y=as.numeric(GBM_train[,108])
cv.fit<-cv.glmnet(data.matrix(x) ,y,family="binomial") 
plot(cv.fit)

fit<-glmnet(data.matrix(x),y,family="binomial", nlambda=500, alpha=1)
# print(fit)
# plot(fit)
plot(fit, xvar="lambda", label=TRUE)



coefficients<-unlist(coef(fit,s=cv.fit$lambda.min))
Active.Index<-which(coefficients!=0) #系数不为0的特征索引
Active.feature<-rownames(coefficients)[which(coefficients!=0)]
Active.coefficients<-coefficients[Active.Index] #系数不为0的特征系数值

  anntrain<-GBM_train[,c(Active.Index-1, ncol(data_MF_exp_1))]
  anntest_R_G<-GBM_test[,c(Active.Index-1, ncol(data_MF_exp_1))]
  LGG_anntest_R_G<-LGG_data_MF_exp_1[,c(Active.Index-1, ncol(LGG_data_MF_exp_1))]
  
  
  ir.nn2<-nnet(group~., data = anntrain,linout=F,size=10,decay=0.001,maxit=100,trace=F, MaxNWts = 100000)
  pre_ann_R_G<-predict(ir.nn2, anntest_R_G, type = "class")
  LGG_pre_ann_R_G<-predict(ir.nn2, LGG_anntest_R_G, type = "class")
  
  
  pre_ann_R_G<-ifelse(pre_ann_R_G=="HIGH", 1, 0)
  anntest_R_G$group<-ifelse(anntest_R_G$group=="HIGH", 1, 0)
  
  LGG_pre_ann_R_G<-ifelse(LGG_pre_ann_R_G=="HIGH", 1, 0)
  LGG_anntest_R_G$group<-ifelse(LGG_anntest_R_G$group=="HIGH", 1, 0)
  
  GBM_auc_R_G <- roc(anntest_R_G$group,pre_ann_R_G)
  LGG_auc_R_G <- roc(LGG_anntest_R_G$group,LGG_pre_ann_R_G)

  GBM<- ggroc(list(GBM_auc_R=GBM_auc_R, GBM_auc_G=GBM_auc_G, GBM_auc_R_G= GBM_auc_R_G), legacy.axes = TRUE)+ 
    theme_bw() + # 更换黑白主题，默认为theme_grey() 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    scale_colour_manual(values=c("#AACF94", "#6A84AA", "#CC6E6E")) + 
    annotate("text", x=0.6, y=0.2, label=paste("AUC of R = ", round(GBM_auc_R$auc, 3), sep = ""),size=3) +
    annotate("text", x=0.6, y=0.3, label=paste("AUC of G = ", round(GBM_auc_G$auc, 3), sep = ""),size=3) +
    annotate("text", x=0.6, y=0.4, label=paste("AUC of R+G = ", round(GBM_auc_R_G$auc, 3), sep = ""),size=3)
  GBM
  LGG <- ggroc(list(LGG_auc_R=LGG_auc_R, LGG_auc_G=LGG_auc_G, LGG_auc_R_G= LGG_auc_R_G), legacy.axes = TRUE)+ 
    theme_bw() + # 更换黑白主题，默认为theme_grey() 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    scale_colour_manual(values=c("#AACF94", "#6A84AA", "#CC6E6E")) + 
    annotate("text", x=0.6, y=0.2, label=paste("AUC of R = ", round(LGG_auc_R$auc, 3), sep = ""),size=3) +
    annotate("text", x=0.6, y=0.3, label=paste("AUC of G = ", round(LGG_auc_G$auc, 3), sep = ""),size=3) +
    annotate("text", x=0.6, y=0.4, label=paste("AUC of R+G = ", round(LGG_auc_R_G$auc, 3), sep = ""),size=3)
  LGG
