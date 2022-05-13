


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
write.table(TIgene_HIGH,"E:\\drug\\TIGENE_HIGHname_21.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


TIgene_HIGH<-read.table("E:\\drug\\TIGENE_HIGHname_21.txt",header=T,sep = "\t", quote = "")
library(clusterProfiler)

eg <- bitr(TIgene_HIGH[,1], 
           fromType="SYMBOL", 
           toType="ENSEMBL",
           OrgDb="org.Hs.eg.db")
##GDSC GBM 细胞系的21个基因的表达
GDSC_geneexp<-read.table("E:\\drug\\GDSC\\sanger1018_brainarray_ensemblgene_rma.txt",header=F,sep = "\t", quote = "")
GDSC_GBM_celllinename<-read.table("E:\\drug\\GDSC\\GBM_cellline_names.txt",header=T,sep = "\t", quote = "")
index_clnames<-match(GDSC_GBM_celllinename[,2],GDSC_geneexp[1,])
GDSC_celllines_geneexp<-GDSC_geneexp[,c(1,index_clnames[-which(is.na(index_clnames))])]
TIGENE_HIGH_celllines_GBM<-merge(eg, GDSC_celllines_geneexp, by.x="ENSEMBL",by.y="V1")
colnames(TIGENE_HIGH_celllines_GBM)<-c("ENSEMBL","SYMBOL",as.character(GDSC_GBM_celllinename[,1][-which(is.na(index_clnames))]))

dim(TIGENE_HIGH_celllines_GBM)
#筛选GBM样本 TIGENE的exp
TIGENE_HIGH_samples_GBM<-merge(as.data.frame(TIgene_HIGH), HIGH_LOW,by.x="x",by.y="gene")[,1:40]
#TIGENE_LOW_samples_GBM<-merge(as.data.frame(TIgene_LOW), HIGH_LOW,by.x="x",by.y="gene")[,c(1,41:89)]
TIGENE_HIGH_samples_GBM<-TIGENE_HIGH_samples_GBM[match(TIGENE_HIGH_celllines_GBM[,2], TIGENE_HIGH_samples_GBM[,1]),]


#标准化
TIGENE_HIGH_samples_GBM<-cbind(TIGENE_HIGH_samples_GBM[,1], scale(TIGENE_HIGH_samples_GBM[,-1],center=T,scale=T))
TIGENE_HIGH_celllines_GBM<-cbind(TIGENE_HIGH_celllines_GBM[,2], scale(TIGENE_HIGH_celllines_GBM[,-(1:2)],center=T,scale=T))



#皮尔森相关系数r，相异性矩阵1-r
samples_celllines_HIGH<-c()
max<-matrix(100,ncol=ncol(TIGENE_HIGH_samples_GBM), nrow=(ncol(TIGENE_HIGH_celllines)-1))
max[1,]<-colnames(TIGENE_HIGH_samples_GBM)
max[,1]<-colnames(TIGENE_HIGH_celllines_GBM)
for(i in 2:ncol(TIGENE_HIGH_samples_GBM)){
	for(j in 2:ncol(TIGENE_HIGH_celllines_GBM)){
	corRESULT<-cor.test(TIGENE_HIGH_samples_GBM[,i],as.numeric(TIGENE_HIGH_celllines_GBM[,j]),method="pearson")
	p_value<-corRESULT[3]
	r_value<-corRESULT[4]
	max[j,i]<-1-as.numeric(r_value)
	if(r_value>0.8 & p_value<0.01){
samples_celllines_HIGH_1<-cbind(colnames(TIGENE_HIGH_samples_GBM)[i],colnames(TIGENE_HIGH_celllines_GBM)[j],r_value)
samples_celllines_HIGH<-rbind(samples_celllines_HIGH,samples_celllines_HIGH_1)
}
	
}
	
}
#相关系>0.8的样本-细胞系对
#a<-as.data.frame(unlist(unique(samples_celllines_HIGH[,2])))

#层次聚类
dt1<-max[-1,-1]
colnames(dt1)<-max[1,-1]
rownames(dt1)<-max[-1,1]
result<- dist(dt1)
result_hc <- hclust(d = result, method = "ward.D")
#dev.new()
plot(result_hc,ylab="distance",xlab="cell lines")

re <- rect.hclust(result_hc, k = 2)
iris.id <- cutree(result_hc, 2)
unGBM_Cls<-as.data.frame(iris.id[which(iris.id==1)])
GBM_Cls<-as.data.frame(iris.id[which(iris.id==2)])
unGBM_Cls<-cbind(rownames(unGBM_Cls), unGBM_Cls)
GBM_Cls<-cbind(rownames(GBM_Cls), GBM_Cls)
colnames(unGBM_Cls)<-c("celllines_names","group")
colnames(GBM_Cls)<-c("celllines_names","group")

write.table(GBM_Cls,"E:\\drug\\CCLE\\gene expression/cluster2_result.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


##不同组细胞系之间  每种药物IC50-Z_score值  做T-test 
CCLE_Drug_IC50<-read.table("E:\\drug\\CCLE/ccle_broad_2019/data_drug_treatment_ZSCORE.txt",header=T,sep = "\t", quote = "")
CCLE_Drug_AUC<-read.table("E:\\drug\\CCLE/ccle_broad_2019/data_drug_treatment_AUC.txt",header=T,sep = "\t", quote = "")

CCLE_Drug_IC50_1<-cbind(colnames(CCLE_Drug_ZSCORE),as.data.frame(t(CCLE_Drug_ZSCORE)))
CCLE_Drug_AUC_1<-cbind(colnames(CCLE_Drug_IC50),as.data.frame(t(CCLE_Drug_IC50)))

colnames(CCLE_Drug_IC50_1)<-c("celllines_names",as.character(CCLE_Drug_ZSCORE[,1]))
colnames(CCLE_Drug_AUC_1)<-c("celllines_names",as.character(CCLE_Drug_IC50[,1]))


unGBM_Drug<-merge(unGBM_Cls,CCLE_Drug_IC50_1,by.x="celllines_names",by.y="celllines_names" )
GBM_Drug<-merge(GBM_Cls,CCLE_Drug_IC50_1,by.x="celllines_names",by.y="celllines_names" )

unGBM_AUC<-merge(unGBM_Cls,CCLE_Drug_AUC_1,by.x="celllines_names",by.y="celllines_names" )
GBM_AUC<-merge(GBM_Cls,CCLE_Drug_AUC_1,by.x="celllines_names",by.y="celllines_names" )


drug_p_out<-c()
for(h in 3:ncol(GBM_Drug)){
	unGBM<-unGBM_Drug[which(is.na(unGBM_Drug[,h])==0),h]
	GBM<-GBM_Drug[which(is.na(GBM_Drug[,h])==0),h]
	unGBM_1<-unGBM_IC50[which(is.na(unGBM_IC50[,h])==0),h]
	GBM_1<-GBM_IC50[which(is.na(GBM_IC50[,h])==0),h]
	
	drug_re_p<-t.test(as.numeric(as.matrix(unGBM)),as.numeric(as.matrix(GBM )))[[3]]
	fc<-mean(as.numeric(as.matrix(GBM )))-mean(as.numeric(as.matrix(unGBM)))
	fc_IC50<-mean(as.numeric(as.matrix(GBM_1 )))-mean(as.numeric(as.matrix(unGBM_1)))
	drugi<-c(colnames(unGBM_Drug)[h],drug_re_p,fc,fc_IC50)
	drug_p_out<-rbind(drug_p_out, drugi)



}

colnames(drug_p_out)<-c("ENTITY_STABLE_ID","p","fc_auc","fc_IC50")
drug_p_fdr<-as.data.frame(cbind(drug_p_out,p.adjust(as.numeric(drug_p_out[,2]),method="fdr")))
drug_p_5<-drug_p_fdr[as.matrix(drug_p_fdr$p) < 0.05 & as.matrix(drug_p_fdr$fc_auc)>0 & as.matrix(drug_p_fdr$fc_IC50)<0,]
drug_p_5<-drug_p_fdr[as.matrix(drug_p_fdr$p) < 0.05,]



unGBM_Drug_IC50<-merge(unGBM_Cls,CCLE_Drug_IC50_1,by.x="celllines_names",by.y="celllines_names")
GBM_Drug_IC50<-merge(GBM_Cls,CCLE_Drug_IC50_1,by.x="celllines_names",by.y="celllines_names" )

CELLLINES_Drug<-merge(drug_p_5,CCLE_Drug_IC50,by.x="ENTITY_STABLE_ID",by.y="ENTITY_STABLE_ID" )




drug<-as.data.frame(table(CELLLINES_Drug[which(CELLLINES_Drug$IC50..uM.<8),4]))


write.table(drug,"E:\\drug\\drug_result.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


