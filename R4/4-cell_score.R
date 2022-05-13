
##################### 	score
#####  A： cnv  PCC
gene_CopyNumber<-read.table("E:\\GBMdata\\GENE EXP\\CNV\\Gistic2_CopyNumber_Gistic2_all_data_by_genes\\Gistic2_CopyNumber_Gistic2_all_data_by_gene.txt",header=T,sep = "\t", quote = "", fill=T)
gene_annotation<-read.table("D:\\part-time job\\all_gene_annotation.txt",header=F,sep = "\t")
gene_annotation<-as.data.frame(gene_annotation[,2])
colnames(gene_annotation)<-"gene"
group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
HIGH<-group_name[group_name$group=="HIGH",1]


HIGH_index_out<-c()
for(i in 1:length(HIGH)){
  HIGH_index1<-grep(unlist(strsplit(as.character(HIGH[i]),"-"))[3],colnames(gene_CopyNumber))
  if(length(HIGH_index1)!=0){
    HIGH_index_out<-c(HIGH_index_out, HIGH_index1)
  }
}
length(HIGH_index_out)
HIGH_CopyNumber<-gene_CopyNumber[,c(1,HIGH_index_out)]
HIGH_CopyNumber1<-HIGH_CopyNumber[,-1]

HIGH_mean<-cbind(as.character(HIGH_CopyNumber[,1]), apply(HIGH_CopyNumber1,1,mean))
colnames(HIGH_mean)<-c("gene", "mean")
HIGH_mean_gene<-merge(gene_annotation, HIGH_mean, by.x="gene", by.y="gene")
head(HIGH_mean_gene)


ccle_log2CNA<-read.table("E:\\drug\\CCLE\\cellline_ccle_broad\\cellline_ccle_broad\\data_log2CNA.txt",header=T,sep = "\t", quote = "", fill=T, check.names=F)
ccle_log2CNA_CNS<-ccle_log2CNA[,c(1,grep("CENTRAL_NERVOUS_SYSTEM",colnames(ccle_log2CNA)))]

HIGH_all_gene<-merge(ccle_log2CNA_CNS, HIGH_mean_gene, by.x="Hugo_Symbol", by.y="gene")

head(HIGH_all_gene)

A<-c()
for(j in 2:(ncol(HIGH_all_gene)-1)){
  r<-unlist(cor(HIGH_all_gene[,j], as.numeric(as.matrix(HIGH_all_gene[,ncol(HIGH_all_gene)]))))
  A_1<-cbind(colnames(HIGH_all_gene)[j],r)
  A<-rbind(A, A_1)
}
A<-as.data.frame(cbind(A[,1], as.numeric(as.matrix(A[,2]))))

####### PCC   density  
plot(density(A$V2, bw = "sj"), main="cnv--PCC")
rug(A$V2)

#####  B： TP53 or PTEN    D： 7 genes 
gene_mc3<-read.table("E:\\GBMdata\\GENE EXP\\somatic mutation\\GBM_mc3_gene_level.txt",header=T,sep = "\t", quote = "", fill=T)
group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
HIGH<-group_name[group_name$group=="HIGH",1]
LOW<-group_name[group_name$group=="LOW",1]

HIGH_index_out<-c()
for(i in 1:length(HIGH)){
  HIGH_index1<-grep(unlist(strsplit(as.character(HIGH[i]),"-"))[3],colnames(gene_mc3))
  if(length(HIGH_index1)!=0){
    HIGH_index_out<-c(HIGH_index_out, HIGH_index1)
  }
}
HIGH_mc3<-gene_mc3[,c(1,HIGH_index_out)]
HIGH_mc31<-HIGH_mc3[,-1]
HIGH_mc3$p<-apply(HIGH_mc31,1,sum)/length(HIGH_index_out)
HIGH_mc3<-HIGH_mc3[order(-HIGH_mc3$p),]
head(HIGH_mc3)
HIGH_gene<-HIGH_mc3[1,1]      #PTEN

LOW_index_out<-c()
for(i in 1:length(LOW)){
  LOW_index1<-grep(unlist(strsplit(as.character(LOW[i]),"-"))[3],colnames(gene_mc3))
  if(length(LOW_index1)!=0){
    LOW_index_out<-c(LOW_index_out, LOW_index1)
  }
}
LOW_mc3<-gene_mc3[,c(1,LOW_index_out)]
LOW_mc31<-LOW_mc3[,-1]
LOW_mc3$p<-apply(LOW_mc31,1,sum)/length(LOW_index_out)
LOW_mc3<-LOW_mc3[order(-LOW_mc3$p),]
head(LOW_mc3)
LOW_gene<-LOW_mc3[1,1]      #PTEN

ccle_mutations<-read.table("E:\\drug\\CCLE\\cellline_ccle_broad\\cellline_ccle_broad\\data_mutations_extended.txt",header=T,sep = "\t", quote = "", fill=T)
ccle_mutations1<-unique(ccle_mutations[grep("CENTRAL_NERVOUS_SYSTEM", ccle_mutations$Tumor_Sample),c(17, 1,10, 11)])
index_1<-which(ccle_mutations1$Variant_Classification=="RNA" |
                 ccle_mutations1$Variant_Classification=="In_Frame_Del" |
                 ccle_mutations1$Variant_Classification=="Silent" |
                 ccle_mutations1$Variant_Classification=="3'UTR" |
                 ccle_mutations1$Variant_Classification=="5'UTR" |
                 ccle_mutations1$Variant_Classification=="Intron" |
                 ccle_mutations1$Variant_Classification=="3'Flank" |
                 ccle_mutations1$Variant_Classification=="5'Flank" |
                 ccle_mutations1$Variant_Classification=="In_Frame_Ins"
)

ccle_mutations2<-ccle_mutations1[-index_1,]
B<-ccle_mutations2[which(ccle_mutations2$Hugo_Symbol=="PTEN" ),]






##############  bar plot
high_gene<-HIGH_mc3[which(HIGH_mc3[,ncol(HIGH_mc3)]>0.1),c(1,ncol(HIGH_mc3))]
low_gene<-LOW_mc3[which(LOW_mc3[,ncol(LOW_mc3)]>0.1),c(1,ncol(LOW_mc3))]
all_gene<-unique(low_gene[,1], high_gene[,1])
ccle_gene_out<-c()
for(i in 1:length(all_gene)){
  number_samples<-length(unique(ccle_mutations2[which(ccle_mutations2[,2]==as.character(all_gene[i]) ),1]))
  if(number_samples!=0){
    ccle_gene<-cbind(as.character(all_gene[i]), as.numeric(number_samples)/length(unique(ccle_mutations2[,1])))
    ccle_gene_out<-rbind(ccle_gene, ccle_gene_out)
  }
  
}
ccle_gene_out1<-as.data.frame(ccle_gene_out)
ccle_high_out1<-merge(ccle_gene_out1, high_gene, by.x = "V1", by.y = "sample")
ccle_low_out1<-merge(ccle_gene_out1, low_gene, by.x = "V1", by.y = "sample")
Celllines<-cbind(ccle_high_out1[,1:2], "Celllines")
colnames(Celllines)<-c("gene","percentage","tissue")
samples<-cbind(as.character(ccle_high_out1[,1]),as.character(ccle_high_out1[,3]), "samples")
colnames(samples)<-c("gene","percentage","tissue")
ccle_high_2<-as.data.frame(rbind(Celllines,samples)) 

Celllines1<-cbind(ccle_low_out1[,1:2], "Celllines")
colnames(Celllines1)<-c("gene","percentage","tissue")
samples1<-cbind(as.character(ccle_low_out1[,1]),as.character(ccle_low_out1[,3]), "samples")
colnames(samples1)<-c("gene","percentage","tissue")
ccle_low_2<-as.data.frame(rbind(Celllines1,samples1)) 

ggplot(as.data.frame(ccle_high_2), aes(x = gene, y = as.numeric(as.matrix(percentage)), fill = tissue)) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+
  scale_fill_manual(values=c("#8FC4AE", "#C5322C")) +labs(x="Genes",y="Percentage",title="HIGH_group",fill="tissue_type")

ggplot(as.data.frame(ccle_low_2), aes(x = gene, y = as.numeric(as.matrix(percentage)), fill = tissue)) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+
  scale_fill_manual(values=c("#8FC4AE", "#5A5AFF")) +labs(x="Genes",y="Percentage",title="LOW_group",fill="tissue_type")


##############  bar plot  2.0
high_gene<-HIGH_mc3[which(HIGH_mc3[,ncol(HIGH_mc3)]>0.1),c(1,ncol(HIGH_mc3))]
low_gene<-LOW_mc3[which(LOW_mc3[,ncol(LOW_mc3)]>0.1),c(1,ncol(LOW_mc3))]
all_gene<-unique(low_gene[,1], high_gene[,1])
ccle_gene_out<-c()
for(i in 1:length(all_gene)){
  number_samples<-length(unique(ccle_mutations2[which(ccle_mutations2[,2]==as.character(all_gene[i]) ),1]))
  if(number_samples!=0){
    ccle_gene<-cbind(as.character(all_gene[i]), as.numeric(number_samples)/length(unique(ccle_mutations2[,1])))
    ccle_gene_out<-rbind(ccle_gene, ccle_gene_out)
  }
  
}
ccle_gene_out1<-as.data.frame(ccle_gene_out)
ccle_high_out1<-merge(ccle_gene_out1, high_gene, by.x = "V1", by.y = "sample", all.x = T)
ccle_high_out1<-merge(ccle_high_out1, low_gene, by.x = "V1", by.y = "sample", all.x = T)
ccle_high_out1[is.na(ccle_high_out1[,3]),3]=0



Celllines1<-cbind(ccle_high_out1[,1:2], "Celllines")
colnames(Celllines1)<-c("gene","percentage","tissue")

samples1<-cbind(as.character(ccle_high_out1[,1]),as.character(ccle_high_out1[,3]), "HIGH")
colnames(samples1)<-c("gene","percentage","tissue")

samples2<-cbind(as.character(ccle_high_out1[,1]),as.character(ccle_high_out1[,4]), "LOW")
colnames(samples2)<-c("gene","percentage","tissue")

ccle_low_2<-rbind(Celllines1,samples1, samples2)


ggplot(as.data.frame(ccle_low_2), aes(x = gene, y = as.numeric(as.matrix(percentage)), fill = tissue)) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+
  scale_fill_manual(values=c("#8FC4AE", "#CB5D4A", "#438AB3")) +labs(x="Genes",y="Percentage",title="",fill="type")+
  theme(panel.background = element_blank())


D<-as.data.frame(ccle_mutations2, HIGH_mc3[which(HIGH_mc3$p>0.1),c(1,ncol(HIGH_mc3))])


#D_gene<-as.data.frame(setdiff(LOW_mc3[which(LOW_mc3$p>0.1),1], HIGH_mc3[which(HIGH_mc3$p>0.1),1]))
#colnames(D_gene)<-"gene"
#D<-merge(ccle_mutations2, D_gene, by.x="Hugo_Symbol", by.y="gene")
##################### D: bar   plot
#HIGHdata<-cbind(merge(HIGH_mc3, D_gene, by.x="sample", by.y="gene")[,c(1,ncol(HIGH_mc3))],"HIGH")
#LOWdata<-cbind(merge(LOW_mc3, D_gene, by.x="sample", by.y="gene")[,c(1,ncol(LOW_mc3))],"LOW")
#colnames(HIGHdata)<-c("gene","percentage","group")
#colnames(LOWdata)<-c("gene","percentage","group")

#D_plot_data<-rbind(HIGHdata,LOWdata)

#ggplot(as.data.frame(D_plot_data), aes(x = gene, y = as.numeric(as.matrix(percentage)), fill = group)) +
# 条形图函数：position设置条形图类型为簇状
# geom_bar(position = "dodge", stat = "identity")+coord_flip()+
# scale_fill_manual(values=c("#CB5D4A", "#438AB3")) +labs(x="Genes",y="Percentage",title="",fill="group")+
# theme(panel.background = element_blank())


#####  C： hypermatution  
setwd("E:\\drug\\CCLE\\CCLE_hybrid_capture1650_hg19_coverage_2012.06.19")
filenames<-dir()
filenames_CNS<-filenames[grep("CENTRAL_NERVOUS_SYSTEM", filenames)]

gene_mutations<-read.table("E:\\drug\\CCLE\\CELL_Mutation_Count_vs_Fraction_of_Genome_Altered.txt",header=T,sep = "\t", quote = "")
a<-gene_mutations[,3]
names(a)<-rownames(gene_mutations)
a1<-matrix(unlist(sapply(a, 
                         myfun<-function(x){
                           INDEX<-grep(x, filenames_CNS)
                           if(length(INDEX)==1){
                             cbind(names(x), INDEX)
                           }	
                         })), ncol=2, byrow = T)

## a1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               :mutation_celllines_names_row  filenames_CNS_row
##  eg:      870                                  23  
##MOGGCCM_CENTRAL_NERVOUS_SYSTEM  MOGGCCM_CENTRAL_NERVOUS_SYSTEM.tumor_depth.wig.txt
MF_cell_out<-c()
for(i in 1:nrow(a1)){
  gene_reads<-read.table(filenames_CNS[as.numeric(a1[i,2])], header=T,sep = "\t", quote = "")
  gene_reads_1<-gene_reads[1:grep("fixedStep chrom=X",gene_reads[,1])[1],]
  gene_reads_2<-gene_reads_1[-grep("fixedStep",gene_reads_1)]
  
  cellline_mutations<-gene_mutations[as.numeric(a1[i,1]),4]
  MF_cell<-cbind(gene_mutations[as.numeric(a1[i,1]),], cellline_mutations/sum(as.numeric(gene_reads_2) >= 1)*1000000)
  MF_cell_out<-rbind(MF_cell_out, MF_cell)
  
  print(i)
  
}
MF_cell_out<-unique(MF_cell_out)
dim(MF_cell_out)

write.table(MF_cell_out,"E:\\drug\\CCLE\\mutation\\gbm_MF_cell_out.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
MF_cell_out<-read.table("E:\\drug\\CCLE\\mutation\\gbm_MF_cell_out.txt",header=T,sep = "\t", quote = "", fill=T)
colnames(MF_cell_out)<-c(colnames(MF_cell_out)[1:5],"MutationFrequency")

ggplot(MF_cell_out, aes(MutationFrequency, CNA.Fraction)) + geom_point()+xlim(0,10)+
  theme(panel.background = element_blank())


##### s


all_data<-c()
s_out<-c()
for(i in 1:nrow(A)){
  B_num<-ifelse(length(which(B[,1]==as.character(A[i,1]) ))==0,0,1)
  D_num<-length(which(D[,2]==as.character(A[i,1])))
  MF<-ifelse(length(which(MF_cell_out[,3]==as.character(A[i,1])))==0,0,MF_cell_out[which(MF_cell_out[,3]==as.character(A[i,1])),6])
  cnv_fre<-ifelse(length(which(MF_cell_out[,3]==as.character(A[i,1])))==0,0,MF_cell_out[which(MF_cell_out[,3]==as.character(A[i,1])),5])
  s<-cbind(unlist(strsplit(as.character(A[i,1]),"_"))[1], as.numeric(as.matrix(A[i,2]))+as.numeric(as.matrix(B_num))-as.numeric(D_num)/6)
  s_out<-rbind(s_out, s)
  all_data1<-cbind(unlist(strsplit(as.character(A[i,1]),"_"))[1], as.numeric(as.matrix(A[i,2])), as.numeric(as.matrix(B_num)), as.numeric(D_num)/6, s[,2], as.numeric(MF) , as.numeric(cnv_fre))
  all_data<-rbind(all_data1, all_data)
  
}

s_out1<-s_out[order(-as.numeric(as.matrix(s_out[,2]) )),]
all_data_1<-all_data[order(-as.numeric(as.matrix(all_data[,5]) )),]
#write.table(s_out1,"E:\\drug\\CCLE\\celllines_score.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

###################################### ccle celllines score bar  plot 
score_data<-as.data.frame(all_data_1)
colnames(score_data)<-c("celllines_name","Corralation","PTEN","specificgenes","score", "MutationFrequency","FGA")
##  1
ggplot(score_data, aes(x = reorder(celllines_name,as.numeric(as.matrix(score))), y = as.numeric(as.matrix(score)))) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+labs(x="celllines_name",y="score")+
  theme(panel.background = element_blank(), panel.border = element_blank())
##  2
ggplot(score_data, aes(x = reorder(celllines_name,as.numeric(as.matrix(score))), y = as.numeric(as.matrix(Corralation)))) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+labs(x="celllines_name",y="Corralation")+
  theme(panel.background = element_blank(), panel.border = element_blank())
##  3
ggplot(score_data, aes(x = reorder(celllines_name,as.numeric(as.matrix(score))), y = as.numeric(as.matrix(PTEN)))) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+labs(x="celllines_name",y="PTEN")+
  theme(panel.background = element_blank(), panel.border = element_blank())+
  geom_text(aes(label=as.numeric(as.matrix(PTEN))),stat="count",vjust=-0.5)
##  4
ggplot(score_data, aes(x = reorder(celllines_name,as.numeric(as.matrix(score))), y = as.numeric(as.matrix(specificgenes)))) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+labs(x="celllines_name",y="specificgenes")+
  theme(panel.background = element_blank(), panel.border = element_blank())
##  5
ggplot(score_data, aes(x = reorder(celllines_name,as.numeric(as.matrix(score))), y = as.numeric(as.matrix(MutationFrequency)))) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+labs(x="celllines_name",y="MutationFrequency")+
  theme(panel.background = element_blank(), panel.border = element_blank())
##  6
ggplot(score_data, aes(x = reorder(celllines_name,as.numeric(as.matrix(score))), y = as.numeric(as.matrix(FGA)))) +
  # 条形图函数：position设置条形图类型为簇状
  geom_bar(position = "dodge", stat = "identity")+coord_flip()+labs(x="celllines_name",y="FGA")+
  theme(panel.background = element_blank(), panel.border = element_blank())

################## case  control celllines

cellline_exp_r0_8<-as.data.frame(sapply(r0_8, my<-function(x){gsub("\\-", "",x)}))
cellline_exp_r8_10<-as.data.frame(sapply(r8_10, my<-function(x){gsub("\\-", "",x)}))

bad_model<-setdiff(cellline_exp_r0_8[,1], s_out1[s_out1[,2] > 1,1])
good_model<-intersect(cellline_exp_r8_10[,1], s_out1[s_out1[,2] > 1,1])


GDSC_DRUG<-read.table("E:\\drug\\GDSC\\GDSC1_DRUG.txt",header=T,sep = "\t", quote = "")
#GDSC_DRUG<-read.table("E:\\drug\\GDSC\\GDSC2_DRUG.txt",header=T,sep = "\t", quote = "")
GDSC_DRUG[,(ncol(GDSC_DRUG)+1)]<-sapply(GDSC_DRUG[,5],  my<-function(x){gsub("\\-", "",x)})
DRUG_NAME<-unique(GDSC_DRUG$DRUG_NAME)
drug_out<-c()
for(i in 1:length(DRUG_NAME)){
  index_drug<-which(GDSC_DRUG$DRUG_NAME==DRUG_NAME[i])
  control<-merge(GDSC_DRUG[index_drug,], as.data.frame(bad_model), by.x="V20", by.y="bad_model")$Z_SCORE
  case<-merge(GDSC_DRUG[index_drug,], as.data.frame(good_model), by.x="V20", by.y="good_model")$Z_SCORE
  p_less<-wilcox.test(case, control, alternative="less",exact=FALSE,correct=FALSE)[3]
  #p_greater<-wilcox.test(case, control, alternative="greater",exact=FALSE,correct=FALSE)[3]
  
  if(p_less < 0.05){
    durg_un<-cbind(as.character(DRUG_NAME[i]), p_less)
    drug_out<-rbind(drug_out, durg_un)
    print(i)
  }
  
}
drug_out

#write.table(drug_out,"E:\\drug\\HIGH_drug_GDSC1.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
#write.table(drug_out,"E:\\drug\\HIGH_drug_GDSC2.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
#################   drug logic50  boxplot
library(ggplot2)
library(ggsignif)
i=239
index_drug<-which(GDSC_DRUG$DRUG_NAME==DRUG_NAME[i])
control<-merge(GDSC_DRUG[index_drug,], as.data.frame(bad_model), by.x="V20", by.y="bad_model")$Z_SCORE
case<-merge(GDSC_DRUG[index_drug,], as.data.frame(good_model), by.x="V20", by.y="good_model")$Z_SCORE
bardata<-as.data.frame(rbind(cbind(case,"case"),cbind(control, "control"))) 
colnames(bardata)<-c("logIC50", "group")
ggplot(bardata,aes(x=group,as.numeric(as.matrix(logIC50)),colour=factor(group)))+
  geom_boxplot(col=c("#CB5D4A","gray"))+
  geom_signif(comparisons = list(c("case","control")),map_signif_level = T,y_position=2.5, test = "wilcox.test", col="black")+
  labs(x="group",y="logIC50", title =as.character(DRUG_NAME[i]))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5)) 


