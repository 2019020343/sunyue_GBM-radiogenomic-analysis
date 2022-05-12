rm(list=ls())
library(survival)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(survminer)
  
setwd("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA")
feature<-read.table("F-S_matrix.txt",sep="\t",head=T)
clidata<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\GBMtest\\PCA\\os_data.txt",sep="\t",head=T)


outdata<-c()
for(i in 2:dim(feature)[2]){
  s1<-unlist(strsplit(unlist(strsplit(colnames(feature)[i],"X",fix = T))[2],".",fix = T)) 
  index0<-grep(s1[1],clidata[,1])
  index2<-grep(s1[2],clidata[,1])
  index1<-intersect(index0,index2)

if(length(index1)!=0){
   index<-max(index1)
   if(clidata$vital_status[index]=="Alive"){
      data<-cbind(as.character( clidata$bcr_patient_barcode[index]),  
			t(feature[,i]), 
			0, 
                  clidata$days_to_last_followup[index])
}else{
      data<-cbind(as.character( clidata$bcr_patient_barcode[index]), 
			t(feature[,i]), 
			1, 
                  clidata$days_to_death[index])

}


outdata<-rbind(outdata,data)
}
}
colnames(outdata)<-c("sample",feature[,1],"status","time")
#write.table(outdata,"os_data1.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

phe<-as.data.frame(outdata)
phe$year=as.numeric(as.matrix(phe$time))/365





## 批量生存分析 使用 logrank test 方法

mySurv=with(phe,Surv(time, status))
log_rank_p <- lapply(2:8, function(i){
 thr=sort(phe$riskout)[round(nrow(phe)*i/10)]
 phe$group=ifelse(phe$riskout> thr,'HIGH','LOW') 
 print(table( phe$group ))
 data.survdiff=survdiff(mySurv~group,data=phe)
 p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
 return(p.val)
}) 
log_rank_p=unlist(log_rank_p)
log_rank_p



i=6
thr=sort(phe$riskout)[round(nrow(phe)*i/10)]
phe$group=ifelse(phe$riskout> thr,'more','less') 
print(table( phe$group ))
sfit <- survfit(Surv(year, status)~group, data=phe) 
ggsurvplot(sfit, conf.int=F, pval=TRUE)


