library(ggplot2)
up<-read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\result\\HT_HG-U133A\\up.txt",header=T,sep = "\t")
down<-read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\result\\HT_HG-U133A\\down.txt",header=T,sep = "\t")
immune1<-read.table("E:\\immune\\immune_gene.txt",header=F,sep = "\t")
up_imm<-as.data.frame(intersect(up[,1], immune1[,1])) 
colnames(up_imm)<-c("v1")
down_imm<-as.data.frame(intersect(down[,1], immune1[,1])) 
colnames(down_imm)<-c("v1")


gene<-read.table("HT_HG-U133A\\HT_HG-U133A_gene_limma_result.txt",header=T,sep = "\t", quote = "")

up_imm1<-merge(UP_imm ,gene, by.x = "v1", by.y = "gene")[1:89]
down_imm1<-merge(down_imm ,gene, by.x = "v1", by.y = "gene")[1:89]

mat<-down_imm1[,2:89]
rownames(mat)<-down_imm1[,1]

pheatmap(mat,scale ="row",
         fontsize = 5,
         fontsize_row = 10, 
         cluster_cols = F,
         color = colorRampPalette(c("blue", "white", "red"))(200), 
         border_color = "NA", 
         show_rownames = T, 
         show_colnames = F, 
         height=5
)

down_imm<-intersect(down[,1], immune1[,1])
up_imm1<-as.data.frame(up_imm[,-1])
rownames(up_imm1)<-up_imm[,1]
colnames(up_imm1)<-c("v1")



length(up_imm)
length(down_imm)



TIGENE<-read.table("E:\\drug\\TIGENE_HIGHname_21.txt",header=T,sep = "\t")
TIGENE_imm<-intersect(TIGENE[,4], immune1[,1])

cg_gene<-read.table("E:\\GBMdata\\GENE EXP\\Methylation\\HumanMethylation450\\cg_gene_match_result.txt",header=T,sep = "\t")
hyper<-cg_gene[cg_gene$diff_moreCless>0.2,1]
hypo<-cg_gene[cg_gene$diff_moreCless< -0.2,1]

length(hyper)
length(hypo)

hyper_imm<-intersect(hyper, immune1[,1])
hypo_imm<-intersect(hypo, immune1[,1])
length(hyper_imm)
length(hypo_imm)

intersect(down_imm, hyper_imm)
intersect(up_imm, hypo_imm)

cnv_gene<-read.table("E:\\GBMdata\\GENE EXP\\CNV\\cnv_out_fdr.txt",header=T,sep = "\t")
cnv_gene_sig<-cnv_gene[cnv_gene$p<0.01,]
cnv_imm<-intersect(cnv_gene_sig[,1], immune1[,1])
length(cnv_imm)

