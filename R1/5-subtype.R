library(reshape2)

#读取数据，按总丰度大小排序，并重排为便于 ggplot2 识别的样式
subtype <- read.delim('E:\\GBMdata\\subtype.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group_name<-read.table("E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA\\group.txt",header=T,sep = "\t")
HIGH<-group_name[which(group_name$group=="HIGH"),c(1,12)]
LOW<-group_name[which(group_name$group=="LOW"),c(1,12)]


colnames(HIGH)<-c("samples", "status")
colnames(LOW)<-c("samples", "status")


HIGHdata<-as.data.frame(merge(HIGH, subtype, by.x = "samples", by.y = "Case",sort = F))
LOWdata<-as.data.frame(merge(LOW, subtype, by.x = "samples", by.y = "Case",sort = F))


HIGH_status<-as.data.frame(table(HIGHdata[,2])) 
LOW_status<-as.data.frame(table(LOWdata[,2])) 
HIGH_IDH<-as.data.frame(table(HIGHdata[,3])) 
LOW_IDH<-as.data.frame(table(LOWdata[,3])) 
HIGH_codeletion<-as.data.frame(table(HIGHdata[,4])) 
LOW_codeletion<-as.data.frame(table(LOWdata[,4])) 
HIGH_MGMT<-as.data.frame(table(HIGHdata[,5])) 
LOW_MGMT<-as.data.frame(table(LOWdata[,5])) 
HIGH_Original_Subtype<-as.data.frame(table(HIGHdata[,6])) 
LOW_Original_Subtype<-as.data.frame(table(LOWdata[,6])) 
HIGH_Grade<-as.data.frame(table(HIGHdata[,7])) 
LOW_Grade<-as.data.frame(table(LOWdata[,7])) 



library(ggplot2)
library(ggpubr)
##绘制饼图时，若展示为空心的“饼环”样式，则也可视为一种圆环图，以 time1 为例作图展示
################  OS  ###############
p1<-ggdonutchart(HIGH_status, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c( "#1D1D1A") 
)
p2<-ggdonutchart(LOW_status, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#CBCBCB", "#1D1D1A") 
)
################  IDH  ################
p3<-ggdonutchart(HIGH_IDH, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#E3DA35","#EAC868") 
)
p4<-ggdonutchart(LOW_IDH, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#E3DA35","#EAC868")
)
################  codeletion  ################
p5<-ggdonutchart(HIGH_codeletion, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#F2B27E")
)
p6<-ggdonutchart(LOW_codeletion, "Freq",
                  label = "Var1",                               
                  fill = "Var1",                            
                  color = "white",                                
                  palette = c("#F2B27E")
)
################  MGMT  ################
p7<-ggdonutchart(HIGH_MGMT, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#AACF94", "#56A02A")
)
p8<-ggdonutchart(LOW_MGMT, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#AACF94", "#56A02A")
)
################  Original_Subtype  ################
p9<-ggdonutchart(HIGH_Original_Subtype, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#D6DCF2", "#B0B7D5", "#8994C0", "#6270AC", "#3B4D97")
)
p10<-ggdonutchart(LOW_Original_Subtype, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#D6DCF2", "#B0B7D5", "#8994C0", "#6270AC", "#3B4D97")
)

################  Grade  ################
p11<-ggdonutchart(HIGH_Grade, "Freq",
                 label = "Var1",                               
                 fill = "Var1",                            
                 color = "white",                                
                 palette = c("#C06261")
)
p12<-ggdonutchart(LOW_Grade, "Freq",
                  label = "Var1",                               
                  fill = "Var1",                            
                  color = "white",                                
                  palette = c("#C06261")
)



################  plot ###############
#install.packages("Rmisc")
library(gridExtra)
#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=6)
grid.arrange(p1,p3,p5,p7,p9,p11,p2,p4,p6,p8,p10,p12,ncol=6)



