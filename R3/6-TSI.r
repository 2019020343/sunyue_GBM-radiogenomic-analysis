

##################TSI
rm(list=ls())
gene_out<-read.table("D:\\part-time job\\1.12\\1-1gene_up5\\gene.txt",header=T,sep = "\t", quote = "\"'")
immcelltype1<-read.table("D:\\part-time job\\1.12\\1-1gene_max\\all.txt",header=T,sep = "\t", quote = "\"'")
max<-apply(immcelltype1[,-(1:2)],1,max)
immcelltype12<-cbind(immcelltype1[,1:2], immcelltype1[,-(1:2)]/max)
immcelltype_19<-cbind(
		apply(immcelltype12[,3:10],1,max),    #1 B.cell
		apply(immcelltype12[,11:13],1,max),  #2 CD34..Hematopoetic.Stem.cells
		apply(immcelltype12[,14:21],1,max),  #3 gamma.delta.T.cell
		apply(immcelltype12[,22:30],1,max),  #4 DC_Control
		apply(immcelltype12[,31:32],1,max),  #5 gamma.delta.T.cell.y
		apply(immcelltype12[,33:38],1,max),  #6 Immature.DCs 
		apply(immcelltype12[,39:42],1,max),  #7 mast.cell
		apply(immcelltype12[,43:46],1,max),  #8 mDC
		apply(immcelltype12[,47:52],1,max),  #9 Monocyte
		apply(immcelltype12[,53:59],1,max),  #10 Neutrophils 
		apply(immcelltype12[,60:67],1,max),  #11 NK_Donor 
		apply(immcelltype12[,68:80],1,max),  #12 pDC_healthy
		apply(immcelltype12[,81:84],1,max),  #13 Resting_CD4.T.cell 
		apply(immcelltype12[,85:88],1,max),  #14 Resting_NK
		apply(immcelltype12[,89:94],1,max),  #15 Resting_NKT
		apply(immcelltype12[,95:97],1,max),  #16 T_cell.CD4.
		apply(immcelltype12[,98:101],1,max), #17 T_cell.CD8._naive
		apply(immcelltype12[,102:107],1,max),#18 T_cell.CD8. 
		apply(immcelltype12[,108:116],1,max)) #19 Th17
immcelltype_19_1<-1-immcelltype_19

a<-as.data.frame(table(gene_out))
gene_1_19<-a[(a$Freq==1 | a$Freq==19),1]
TSI_data<-c()


for(j in 1:length(gene_1_19)){
	index<-which(immcelltype12[,1]==as.character(gene_1_19[j]))
	if(length(index)==1){
	TSI_j<-cumsum(immcelltype_19_1[index,])[19]/18
	TSI_data1<-cbind(immcelltype12[index,],TSI_j)
	TSI_data<-rbind(TSI_data,TSI_data1)
}else{
print(gene_1_19[j])

}

	
}

write.table(TSI_data,"D:\\part-time job\\1.12\\2MGN\\TSI18.txt" ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


hkGENE<-TSI_data[TSI_data$TSI_j<0.2,]

write.table(hkGENE,"D:\\part-time job\\1.12\\2MGN\\hkGENE18.txt"  ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
icsgene<-TSI_data[TSI_data$TSI_j>0.5,]
write.table(icsgene,"D:\\part-time job\\1.12\\2MGN\\icsgene18.txt"  ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)



MGN<-read.table("D:\\part-time job\\1.12\\2MGN\\limma_result.txt",header=T,sep = "\t", quote = "\"'")
MGN<-MGN[MGN$adj.P.Val<0.01,]
up<-unique(merge(MGN[as.numeric(MGN$logFC)>1,],hkGENE,by="gene"))
dim(up)

write.table(up,"D:\\part-time job\\1.12\\2MGN\\hkGENE_upMGN18.txt"  ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)

down<-unique(merge(MGN[as.numeric(MGN$logFC)<(-1),],hkGENE,by="gene"))
dim(down)


write.table(down,"D:\\part-time job\\1.12\\2MGN\\hkGENE_dowNMGN18.txt"  ,col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


