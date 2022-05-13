rm(list=ls())
#########################   GBM celllines
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

## a1:mutation_celllines_names_row  filenames_CNS_row
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

plot(MF_cell_out[,6], MF_cell_out[,5])
#########################   ov   celllines

rm(list=ls())
##    celllines
setwd("E:\\drug\\CCLE\\CCLE_hybrid_capture1650_hg19_coverage_2012.06.19")
filenames<-dir()
filenames_CNS<-filenames[grep("OVARY", filenames)]

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

## a1:mutation_celllines_names_row  filenames_CNS_row
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
dim(MF_cell_out)
write.table(MF_cell_out,"E:\\drug\\CCLE\\mutation\\OV_MF_cell_out.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)





###     OV  tcga
rm(list=ls())
mc3<-read.table("E:\\drug\\tcgaov\\OV_mc3.txt",header=T,sep = "\t", quote = "", fill=T)
gene_mc3<-read.table("E:\\drug\\tcgaov\\OV_mc3_gene_level.txt",header=T,sep = "\t", quote = "", fill=T)
gene_position<-read.table("E:\\drug\\gene_position.txt",header=F,sep = "\t", quote = "", fill=T)
gene_mutations<-read.table("E:\\drug\\tcgaov\\Mutation_Count_vs_Fraction_of_Genome_Altered (1).txt",header=T,sep = "\t", quote = "")

index_1<-which(mc3$effect=="RNA" |
		 mc3$effect=="In_Frame_Del" |
		 mc3$effect=="Silent" |
		 mc3$effect=="3'UTR" |
		 mc3$effect=="5'UTR" |
		 mc3$effect=="Intron" |
		 mc3$effect=="3'Flank" |
		 mc3$effect=="5'Flank" |
		 mc3$effect=="In_Frame_Ins"
)
mc3_1<-mc3[-index_1,]
mc3_1[1,]
sample_name<-substr(colnames(gene_mc3),9,12)

MF_TCGA_OUT<-c()
for(i in 2:length(sample_name)){
	#counts<-gene_mutations[grep(sample_name[i], gene_mutations[,2]),4]
	counts<-length(grep(sample_name[i], mc3_1[,1]))
	if(length(counts)==1){
	name1<-as.data.frame(mc3_1[grep(sample_name[i], mc3_1[,1]),7])
	colnames(name1)<-"gene"
	position1<-merge(name1,gene_position, by.x="gene", by.y="V4" )
	#MF_1<-cbind(gene_mutations[grep(sample_name[i], gene_mutations[,2]),], counts/sum(position1$V3-position1$V2)*1000000)
	MF_1<-cbind(as.character(unique(mc3_1[grep(sample_name[i], mc3_1[,1]),1])), counts, counts/sum(position1$V3-position1$V2)*1000000)
	MF_TCGA_OUT<-rbind(MF_TCGA_OUT, MF_1)
print(i)
}
}


###     GBM  tcga
rm(list=ls())
mc3<-read.table("E:\\GBMdata\\GENE EXP\\somatic mutation\\GBM_mc3.txt",header=T,sep = "\t", quote = "", fill=T)
gene_mc3<-read.table("E:\\GBMdata\\GENE EXP\\somatic mutation\\GBM_mc3_gene_level.txt",header=T,sep = "\t", quote = "", fill=T)
gene_position<-read.table("E:\\drug\\gene_position.txt",header=F,sep = "\t", quote = "", fill=T)
gene_mutations<-read.table("E:\\GBMdata\\GENE EXP\\somatic mutation\\Mutation_Count_vs_Fraction_of_Genome_Altered (1).txt",header=T,sep = "\t", quote = "")

index_1<-which(mc3$effect=="RNA" |
		 mc3$effect=="In_Frame_Del" |
		 mc3$effect=="Silent" |
		 mc3$effect=="3'UTR" |
		 mc3$effect=="5'UTR" |
		 mc3$effect=="Intron" |
		 mc3$effect=="3'Flank" |
		 mc3$effect=="5'Flank" |
		 mc3$effect=="In_Frame_Ins"
)
mc3_1<-mc3[-index_1,]
mc3_1[1,]
sample_name<-substr(colnames(gene_mc3),9,12)

MF_TCGA_OUT<-c()
for(i in 2:length(sample_name)){
	counts<-gene_mutations[grep(sample_name[i], gene_mutations[,2]),4]
	#counts<-length(grep(sample_name[i], mc3_1[,1]))
	if(counts!=0){
	#name1<-as.data.frame(mc3_1[grep(sample_name[i], mc3_1[,1]),7])
	name1<-as.data.frame(gene_mc3[which(gene_mc3[,i]==1),1])
	colnames(name1)<-"gene"
	position1<-merge(name1,gene_position, by.x="gene", by.y="V4" )
	MF_1<-cbind(gene_mutations[grep(sample_name[i], gene_mutations[,2]),], counts/sum(position1$V3-position1$V2)*1000000)
	#MF_1<-cbind(as.character(unique(mc3_1[grep(sample_name[i], mc3_1[,1]),1])), counts, counts/sum(position1$V3-position1$V2)*1000000)
	MF_TCGA_OUT<-rbind(MF_TCGA_OUT, MF_1)
print(i)
}
}


write.table(MF_TCGA_OUT,"E:\\drug\\CCLE\\mutation\\gbm_MF_TCGA_out_1.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)


MF_TCGA_OUT<-read.table("E:\\drug\\CCLE\\mutation\\gbm_MF_TCGA_out.txt",header=T,sep = "\t", quote = "", fill=T)







