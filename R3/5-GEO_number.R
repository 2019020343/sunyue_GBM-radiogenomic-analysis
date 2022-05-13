library(ggplot2) #加载ggplot2包
library(dplyr) #加载dplyr包
##install.packages("ggstatsplot")
library(ggstatsplot) #加载ggstatsplot包
library(RColorBrewer)

#使用ggstatsplot的ggbarstats函数
setwd("E:\\immune")
GSEDATA<-read.table("GSEDATA.txt",header=T,sep = "\t", stringsAsFactors=FALSE)



colnames(GSEDATA)<-c("GSE", "celltype", "samples_number")




GSEDATA%>% 
  ggplot(aes(x =reorder(GSE, -as.numeric(samples_number) )  , y=as.numeric(samples_number), fill = celltype )) + #x轴的分类为clarity，填充颜色为color（J和H）
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#B3D5D7","#65ACC0","#9DC77C","#DBDA7A","#C77A25",
                               "#D3B494","#C2ADCB","#69589C", "#9AC5AB","#90B7E2",
                               "#B96765","#B28AB6","#C58D98", "#9E1A31","#D6C591",
                               "#CEA121","#64A1DF","#B0D0DD","#7399A8")) +
  theme_classic() + #设置主题
  labs(x = "GEO", y = 'samples_number')+  #设置y轴名为‘Percent’
  coord_flip() #旋转坐标轴


