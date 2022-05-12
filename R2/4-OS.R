library(survival)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(survminer)  
setwd("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA")
clusterdata<-read.table("PCA-10.txt",sep="\t",head=T)
clidata<-read.table("GBMLGG_survival.txt",sep="\t",head=T)
clusterdata1<-clusterdata
outdata<-c()
for(i in 1:dim(clusterdata1)[1]){
  s1<-clusterdata1[i,1]
  index1<-grep(s1,clidata[,1])
if(length(index1)!=0){
   index<-max(index1)
   data<-cbind(as.character( clidata$X_PATIENT[index]),  
			clusterdata1[i,-1], 
                  clidata$OS[index], 
			clidata$OS.time[index])

colnames(data)<-c(colnames(clusterdata1) ,"status","time")
outdata<-rbind(outdata,data)
}
}
write.table(outdata,"E:/LGG_kpsos.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


# outdata<-read.table("E:/kpsos.txt",sep="\t",head=T)

cox<-coxph(Surv(as.numeric(as.matrix(time)),status)~pca_1+pca_2+pca_3+pca_4+pca_5+pca_6+pca_7+pca_8+pca_9+pca_10, data = outdata)
cox
ggforest(cox, data = outdata)


m<-as.data.frame(summary(cox)$coef)
b_value<-m[,1]

##PCA 1-7+10 显著
riskdata1<-outdata
riskout<-c()
for(j in 1:dim(riskdata1)[1]){

  risk<-riskdata1[j,2]*b_value[1]+
		riskdata1[j,3]*b_value[2]+
		riskdata1[j,4]*b_value[3]+
		riskdata1[j,5]*b_value[4]+
		riskdata1[j,6]*b_value[5]+
		riskdata1[j,7]*b_value[6]+
		riskdata1[j,8]*b_value[7]+
		riskdata1[j,9]*b_value[8]+
		riskdata1[j,10]*b_value[9]+
		riskdata1[j,11]*b_value[10]
  riskout<-c(riskout,risk)


}
riskoutdata<-cbind(riskdata1,riskout)


cox1<-coxph(Surv(as.numeric(as.matrix(time)),status)~riskout, data =  riskoutdata)
cox1
ggforest(cox1, data = riskoutdata)

# phe<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\risk_group\\group.txt",sep="\t",head=T)

phe<-riskoutdata
phe$year=phe$time/365

## 批量生存分析 使用 logrank test 方法
mySurv=with(phe,Surv(year, status))
log_rank_p <- lapply(2:8, function(i){
 thr=sort(phe$riskout)[round(nrow(phe)*i/10)]
 phe$group=ifelse(phe$riskout> thr,'more','less') 
 print(table( phe$group ))
 data.survdiff=survdiff(mySurv~group,data=phe)
 p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
 return(p.val)
}) 
log_rank_p=unlist(log_rank_p)
log_rank_p



i=8
thr=sort(phe$riskout)[round(nrow(phe)*i/10)]
phe$group=ifelse(phe$riskout> thr,'more','less') 
print(table( phe$group ))
sfit <- survfit(Surv(year, status)~group, data=phe) 
ggsurvplot(sfit, conf.int=F, pval=TRUE)
write.table(phe,"E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\risk_group\\group.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)




fit1 <- survfit(Surv(as.numeric(phe[,13]),as.numeric(phe[,12])) ~phe[,16] )
survdiff(Surv(as.numeric(phe[,13]),as.numeric(phe[,12])) ~phe[,16])

ggsurvplot(fit1,data = riskoutdata,

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

