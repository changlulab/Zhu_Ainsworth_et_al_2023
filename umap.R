rm(list = ls())
options(stringsAsFactors = F)
library(umap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scatterplot3d)
options(rgl.useNULL=TRUE)
library(rgl)
rgl::setupKnitr(autoprint = TRUE)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/glia/")
df = read.csv('All_SCZ_Control_glia.csv')
new_df = df[which(df$padj < 0.05 & abs(df$log2FoldChange) >0.26 ),]
tmp_data = new_df[,8:123]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/cofactors/")
species = read.csv('information_3d_plot.csv')
colnames(species)[1] = 'ID'
colnames(tmp_data) = species$ID
data = t(tmp_data)
dim(data)

type='neuron'
mark='Enhancer'
nn = 12
umap_results <- umap(data,n_neighbors=nn,min_dist=0.1,n_components = 3)
embedding = data.frame(umap_results$layout)
embedding$condition =as.factor(species$condition)
embedding$gender = as.factor(species$gender)
embedding$age = as.factor(species$age)
embedding$PMI = as.factor(species$PMI)
colnames(embedding) = c('V1','V2','V3','condition','gender','age','PMI')
plot3d(embedding[1:116,1:3],type = 's',col = as.integer(embedding$condition),size = 2)

pdf("~/Desktop/test.pdf")
ggplot(embedding,mapping = aes(x=V1,y=V2,col=condition,size=age,label=age,shape = factor(gender))) +
  geom_point(aes(colour = condition)) +
  geom_text(data=subset(embedding,
                        condition == 'untrt_SCZ' | condition == 'trt_SCZ'),
            aes(label=age),
            hjust=0,vjust =0,
            nudge_x = 0.05,
            nudge_y = 0.05)

setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/for reviewers/rebuttal_PointToPoint/clustering/UMAP/DEGs/")

pdf(paste0(mark,'_',type,'_UMAP_full_nn',nn,'.pdf'), height = 8,width = 12)
ggplot(embedding,mapping = aes(x=V1,y=V2,col=condition,size=age,label=age,shape = factor(gender))) +
  geom_point(aes(colour = condition)) +
  geom_text(data=subset(embedding,
                        condition == 'untrt_SCZ' | condition == 'trt_SCZ'),
            aes(label=age),
            hjust=0,vjust =0,
            nudge_x = 0.05,
            nudge_y = 0.05)
dev.off()
write.csv(embedding,paste0(mark,'_',type,'_UMAP_full_nn',nn,'.csv'))
