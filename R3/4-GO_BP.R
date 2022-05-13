#BiocManager::install("topGO")
rm(list = ls())
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(ggplot2)

gene=read.table("E:\\GBMdata\\GENE EXP\\RNA-seq\\result\\HT_HG-U133A\\HT_HG-U133A_gene_limma_result.txt",header=T)
up<-gene[gene$P.Value<0.01   & gene$logFC > 0 ,c(1,90)]
down<-gene[gene$P.Value<0.01   & gene$logFC < 0 ,c(1,90)]
up1<-up[,2]
down1<-down[,2]
names(up1)<-up[,1]
names(down1)<-down[,1]


upsort<-sort(up1,decreasing=T)
downsort<-sort(down1,decreasing=T)


##################  GO

upegseGO <- gseGO(upsort, OrgDb=org.Hs.eg.db,
                  ont='ALL',keyType="SYMBOL",
                  nPerm=1000, minGSSize=10, maxGSSize=500,
                  pvalueCutoff=0.05, verbose=FALSE, by="fgsea")
downegseGO <- gseGO(downsort, OrgDb=org.Hs.eg.db,
                ont='ALL',keyType="SYMBOL",
                nPerm=1000, minGSSize=10, maxGSSize=500,
                pvalueCutoff=0.05, verbose=FALSE, by="fgsea")
upegseGO1<-as.data.frame(upegseGO)
downegseGO1<-as.data.frame(downegseGO)

#################  kegg
gene.tx <- bitr(names(upsort),fromType="SYMBOL",toType=c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
colnames(gene.tx)[1] <- "gene"
gene.tx <- merge(gene.tx,up,by="gene")
FCgenelist <- up$logFC #numeric vector
names(FCgenelist) <- as.character(gene.tx$ENTREZID) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order

upegseKEGG <- gseKEGG(FCgenelist, organism = "hsa",
                      keyType = "kegg",
                  nPerm=100, minGSSize=5, maxGSSize=100,
                  pvalueCutoff=1, pAdjustMethod = "BH", by="fgsea")
upegseKEGG1<-as.data.frame(upegseKEGG)


gene.tx <- bitr(names(downsort),fromType="SYMBOL",toType=c("ENTREZID"),OrgDb = org.Hs.eg.db)
colnames(gene.tx)[1] <- "gene"
gene.tx <- merge(gene.tx,down,by="gene")
FCgenelist <- down$logFC #numeric vector
names(FCgenelist) <- as.character(gene.tx$ENTREZID) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order

downegseKEGG <- gseKEGG(FCgenelist, organism = "hsa",
                      keyType = "ncbi-geneid",
                      nPerm=100, minGSSize=5, maxGSSize=100,
                      pvalueCutoff=1, pAdjustMethod = "BH", by="fgsea")
downegseKEGG1<-as.data.frame(downegseKEGG)


#install.packages("ggupset")
library(ggupset)
upsetplot(downegseGO)
heatplot(downegseGO, foldChange=downsort)


###############################  enrich  go
up_ENTREZID=bitr(up[,1], 'SYMBOL', "ENTREZID", "org.Hs.eg.db")  #将Gene名字转化为基因ID，若输入的是基因SYMBOL需要用该命令转化成基因ID，若输入的是基因ID，则忽略该命令。
down_ENTREZID=bitr(down[,1], 'SYMBOL', "ENTREZID", "org.Hs.eg.db") #将Gene名字转化为基因ID，若输入的是基因SYMBOL需要用该命令转化成基因ID，若输入的是基因ID，则忽略该命令。
up_ENTREZID_1<-merge(up_ENTREZID, up, by.x="SYMBOL", by.y="gene" )
down_ENTREZID_1<-merge(down_ENTREZID, down, by.x="SYMBOL", by.y="gene" )


setwd("E:\\GBMdata\\GENE EXP\\RNA-seq\\result\\HT_HG-U133A\\clusterProfiler")
#################  GO   BP
up_GO_BP<- enrichGO(gene = up_ENTREZID_1[,2],
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)       #2


down_GO_BP<- enrichGO(gene = down_ENTREZID_1[,2],
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)       #27


#data<-rbind(as.data.frame(up_GO_BP), as.data.frame(down_GO_BP))
#write.table(data,"go_bp.txt", sep = "\t", quote = F, row.names = F)
#data<-read.table("go_bp.txt", header = T,sep = "\t", quote = "")

  
upsetplot(down_GO_BP)

genelist<- -down_ENTREZID_1[,3]
names(genelist)<-down_ENTREZID_1[,1]
heatplot(down_GO_BP, showCategory = nrow(down_GO_BP), foldChange=genelist)

upsetplot(up_GO_BP)


upgenelist<- up_ENTREZID_1[,3]
names(upgenelist)<-up_ENTREZID_1[,1]
heatplot(up_GO_BP,  showCategory = nrow(up_GO_BP), foldChange=upgenelist)




data[,10]<-rep(c("up","down"),c(nrow(up_GO_BP), nrow(down_GO_BP)))


ggplot(data, aes(V10, Description)) + 
  geom_tile(aes(fill = p.adjust),colour = "white") + 
  scale_fill_gradient2(limits = c(min(data$p.adjust),max(data$p.adjust)),low = "white", mid="black", high = "red", midpoint = 0)+
  theme_classic()+
  coord_fixed()




#################  GO   MF
up_GO_MF<- enrichGO(gene = up_ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)       #7


down_GO_MF<- enrichGO(gene = down_ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)       #2



#################  GO   cc
up_GO_CC<- enrichGO(gene = up_ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)       #0 

down_GO_CC<- enrichGO(gene = down_ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)       #32