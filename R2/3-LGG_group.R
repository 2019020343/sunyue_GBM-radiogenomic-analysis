################# cox b #################
rm(list=ls())
library("survival")
library("survminer")
setwd("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA")
PCA<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA\\PCA-10.txt",sep="\t",header=T)
b_value<-read.table("cox_b_value.txt",sep="\t",header=T)
a<-as.matrix(PCA[,-1]) 
rownames(a)<-PCA[,1]
risk_score<-a %*% as.matrix(b_value) 
clidata<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA\\GBMLGG_survival.txt",sep="\t",head=T)
risk_score<-cbind(risk_score, "OS", "OS.time")
for (i in 1:nrow(risk_score)) {
  index<-grep(rownames(risk_score)[i], clidata[,1])
  if(length(index)!=0){
    risk_score[i,2]<-clidata$OS[index]
    risk_score[i,3]<-clidata$OS.time[index]
    
  }else{
    index_un<-i
  }

  
}
risk_score<-risk_score[-index_un,]
label<-ifelse(risk_score[,1]>0,"HIGH", "LOW")
risk_score1<-cbind(risk_score, label)
colnames(risk_score1)<-c("risk_score", "status", "time", "group")
phe<-as.data.frame(risk_score1)
phe$year=as.numeric(as.matrix(phe$time))/365
sfit <- survfit(Surv(as.numeric(as.matrix(year)) , as.numeric(status))~group, data=phe) 
ggsurvplot(sfit, conf.int=F, pval=TRUE)
################### weight ###############
rm(list=ls())
library("survival")
library("survminer")
setwd("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA")
feature_os_group<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA\\feature_os_group.txt",sep="\t",header=T)
weight<-read.table("weight.txt",sep="\t",header=T)
b_value<-read.table("cox_b_value.txt",sep="\t",header=T)
feature<-as.matrix(feature_os_group[,2:87]) 
feature<-scale(feature)
rownames(feature)<-feature_os_group[,1]
weight1<-matrix(as.numeric( t(weight[,-1])), ncol = 10, nrow = 86)
a<-feature %*% weight1
clidata<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA\\GBMLGG_survival.txt",sep="\t",head=T)
a<-cbind(a, "OS", "OS.time")
for (i in 1:nrow(a)) {
  index<-grep(rownames(a)[i], clidata[,1])
  if(length(index)!=0){
    a[i,11]<-clidata$OS[index]
    a[i,12]<-clidata$OS.time[index]
  }else{
    index_un<-i
  }
}


a1<-matrix(as.numeric(as.matrix(a)), ncol = 12, nrow = 43)
a1<-as.data.frame(a1)
rownames(a1)<-rownames(a)
colnames(a1)<-c("pca_1","pca_2","pca_3","pca_4","pca_5","pca_6","pca_7","pca_8","pca_9","pca_10", "status", "time")
cox<-coxph(Surv(as.numeric(as.matrix(time)),as.numeric(as.matrix(status)))~pca_1+pca_2+pca_3+pca_4+pca_5+pca_6+pca_7+pca_8+pca_9+pca_10, data =a1)
ggforest(cox, data = a1)

m<-as.data.frame(summary(cox)$coef)
b_value<-m[,1]
riskdata1<-a1
riskout<-c()
for(j in 1:dim(riskdata1)[1]){
  
  risk<-riskdata1[j,1]*b_value[1]+
    riskdata1[j,2]*b_value[2]+
    riskdata1[j,3]*b_value[3]+
    riskdata1[j,4]*b_value[4]+
    riskdata1[j,5]*b_value[5]+
    riskdata1[j,6]*b_value[6]+
    riskdata1[j,7]*b_value[7]+
    riskdata1[j,8]*b_value[8]+
    riskdata1[j,9]*b_value[9]+
    riskdata1[j,10]*b_value[10]
  riskout<-c(riskout,risk)
  
  
}
riskoutdata<-cbind(riskdata1,riskout)
cox1<-coxph(Surv(as.numeric(as.matrix(time)),status)~riskout, data =  riskoutdata)
ggforest(cox1, data = riskoutdata)



label<-ifelse(riskoutdata[,13]>0,"HIGH", "LOW")
risk_score1<-cbind(riskoutdata, label)
colnames(risk_score1)<-c(colnames(riskoutdata),"group")
phe<-as.data.frame(risk_score1)
phe$year=as.numeric(as.matrix(phe$time))/365
sfit <- survfit(Surv(as.numeric(as.matrix(year)) , as.numeric(status))~group, data=phe) 
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,data = phe,
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           # submain = coxdata$cerna[1], 
           legend = "top", 
           legend.title = "risk_score"
           # palette = c("#E7B800", "#2E9FDF","red")
)


feature_label<-cbind(feature, risk_score1$group)##标准化的特征和分组矩阵
write.table(risk_score1,"E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA\\group.txt",sep="\t",quote=FALSE, row.name=T, col.name=T)




