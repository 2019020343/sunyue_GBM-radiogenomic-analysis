rm(list=ls())
library(nnet)
library(glmnet)
library(pROC)
library(plyr)
library(glmnet)
library(ggplot2)

GBM_data_MF<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\0.73\\F_M_matrix_BZ.txt",sep="\t",header=T, row.names = 1)
LGG_data_MF<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\lasso\\0.71\\F_M_matrix_BZ.txt",sep="\t",header=T, row.names = 1)
out_08<-c()
out_8<-c()
seedset<-sample(1000:2500000, 500, replace = F)
for (j in 1:500) {
set.seed(seedset[j])
n<-length(GBM_data_MF[,1])
samp<-sample(1:n,n/5)
GBM_train<-GBM_data_MF[-samp, ]
GBM_test<-GBM_data_MF[samp, ]

x=GBM_train[,1:86]
x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], function(m) as.numeric(as.character(m)))

y=as.numeric(GBM_train[,87])
cv.fit<-cv.glmnet(data.matrix(x),y,family="binomial") 
#plot(cv.fit)

fit<-glmnet(as.matrix(x),y,family="binomial", nlambda=500, alpha=1)
# print(fit)
# plot(fit)
#plot(fit, xvar="lambda", label=TRUE)
coefficients<-unlist(coef(fit,s=cv.fit$lambda.min))
Active.Index<-which(coefficients!=0) #系数不为0的特征索引
Active.feature<-rownames(coefficients)[which(coefficients!=0)]
Active.coefficients<-coefficients[Active.Index] #系数不为0的特征系数值
if(length(Active.feature)>8){
  anntrain<-GBM_train[,c(Active.Index-1, ncol(GBM_data_MF))]
  anntest<-GBM_test[,c(Active.Index-1, ncol(GBM_data_MF))]
  LGG_anntest<-LGG_data_MF[,c(Active.Index-1, ncol(LGG_data_MF))]
  
  
  ir.nn2<-nnet(group~., data = anntrain,linout=F,size=5,decay=5e-4,maxit=100,trace=F,rang=0.7)
  pre_ann<-predict(ir.nn2, anntest, type = "class")
  LGG_pre_ann<-predict(ir.nn2, LGG_anntest, type = "class")
  
  
  pre_ann<-ifelse(pre_ann=="HIGH", 1, 2)
  anntest$group<-ifelse(anntest$group=="HIGH", 1, 2)
  
  LGG_pre_ann<-ifelse(LGG_pre_ann=="HIGH", 1, 2)
  LGG_anntest$group<-ifelse(LGG_anntest$group=="HIGH", 1, 2)
  
  GBM_auc <- roc(anntest$group,pre_ann )$auc
  LGG_auc <- roc(LGG_anntest$group,LGG_pre_ann)$auc
  
  out_1<-cbind(j,seedset[j],GBM_auc, LGG_auc)
  out_08<-rbind(out_08, out_1)
  print(j)
  
}

}


################################# final  model  ###############
library(nnet)
library(glmnet)
library(pROC)
library(plyr)
library(glmnet)
library(ggplot2)

GBM_data_MF<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\lasso\\0.73\\F_M_matrix_BZ.txt",sep="\t",header=T, row.names = 1)
LGG_data_MF<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\lasso\\0.71\\F_M_matrix_BZ.txt",sep="\t",header=T, row.names = 1)

set.seed(73902)
n<-length(GBM_data_MF[,1])
samp<-sample(1:n,n/5)
GBM_train<-GBM_data_MF[-samp, ]
GBM_test<-GBM_data_MF[samp, ]

x=as.matrix(GBM_train[,1:86])
y=GBM_train[,87]
cv.fit<-cv.glmnet(data.matrix(x),y,family="binomial") 
#plot(cv.fit)

fit<-glmnet(as.matrix(x),y,family="binomial", nlambda=500, alpha=1)
# print(fit)
#plot(fit)
#plot(fit, xvar="lambda", label=TRUE)
coefficients<-unlist(coef(fit,s=cv.fit$lambda.min))
Active.Index<-which(coefficients!=0) #系数不为0的特征索引
Active.feature<-rownames(coefficients)[which(coefficients!=0)]
Active.coefficients<-coefficients[Active.Index] #系数不为0的特征系数值

  anntrain<-GBM_train[,c(Active.Index-1, ncol(GBM_data_MF))]
  anntest_R<-GBM_test[,c(Active.Index-1, ncol(GBM_data_MF))]
  LGG_anntest_R<-LGG_data_MF[,c(Active.Index-1, ncol(LGG_data_MF))]
  
  
  ir.nn2<-nnet(group~., data = anntrain,linout=F,size=5,decay=5e-4,maxit=100,trace=F,rang=0.7)
  summary(ir.nn2)
  pre_ann_R<-predict(ir.nn2, anntest_R, type = "class")
  LGG_pre_ann_R<-predict(ir.nn2, LGG_anntest_R, type = "class")
  
  
  pre_ann_R<-ifelse(pre_ann_R=="HIGH", 1, 2)
  anntest_R$group<-ifelse(anntest_R$group=="HIGH", 1, 2)
  
  LGG_pre_ann_R<-ifelse(LGG_pre_ann_R=="HIGH", 1, 2)
  LGG_anntest_R$group<-ifelse(LGG_anntest_R$group=="HIGH", 1, 2)
  
  GBM_auc_R <- roc(anntest_R$group,pre_ann_R )#0.758
  LGG_auc_R <- roc(LGG_anntest_R$group,LGG_pre_ann_R)#0.697
  
  
  
  