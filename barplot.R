#############pre#######
rm(list = ls())
options(stringsAsFactors = F)
library(venneuler)
library(ggplot2)
library(ggpubr)
library(viridis)
library(hrbrthemes)
##########################################MIA_birth_mom#######################################
specie <- c(rep("enh_neun+",2),rep("enh_neun-",2),
            rep("pro_neun+",2),rep("pro_neun-",2),
            rep("DEGs_neun+",2),rep("DEGs_neun-",2))
condition <- rep(c("increase" , "decrease"), 6)
value <- c('1802','1856','1075','1576',
           '14','22','281','494',
           '605','668','298','478')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)
data$order = 1:12
theme <- theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA,linetype = 1, color="black",size=1))

# fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie)) + 
#   geom_bar(position="stack", stat="identity") +
#   labs(x=NULL,y=NULL)+    
#   coord_cartesian(ylim = c(0,1000)) +  
#   theme
# 
# fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie)) + 
#   geom_bar(position="stack", stat="identity") +
#   labs(x=NULL,y=NULL) +  
#   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + 
#   coord_cartesian(ylim = c(1000,3000)) +  
#   theme
# 
# ggarrange(fig2,fig1,heights=c(1/2, 1/2),ncol = 1, nrow = 2,
#           common.legend = T,legend="right",align = "v")

setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/Cell reports_review/figures for revision/barplots/new")
pdf('trt_SCZs vs. control_2.pdf', width = 6, height = 5)
ggplot(data,aes(fill=condition, y=value, x=reorder(specie,order))) +
 geom_bar(position="stack", stat="identity") +
  labs(x=NULL,y=NULL)+
  coord_cartesian(ylim = c(0,4000)) +
  theme
dev.off()
