###################neurons#############################
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/")
factors = read.csv('metadata_RNA_neuron.csv')
row.names(factors) = factors[,1]
factors = factors[,-1]
factors$total.reads = round((factors$total.reads/1000000),2)
untrt_Control1_factors = bind_rows(factors[1:30,],factors[59:88,])
trt_Control2_factors = bind_rows(factors[31:58,],factors[89:116,])
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/neuron/")
SCZ_1 = read.table('SCZ-1-8-neuron.txt',header = T)
SCZ_2 = read.table('SCZ-9-16-neuron.txt',header = T)
SCZ_3 = read.table('SCZ-17-24-neuron.txt',header = T)
SCZ_4= read.table('SCZ-25-29-neuron.txt',header = T)
SCZ_count = data.frame(SCZ_1[,1],SCZ_1[,7:ncol(SCZ_1)],SCZ_2[,21:22],SCZ_2[,7:20],
                       SCZ_3[,7:ncol(SCZ_3)],SCZ_4[,7:ncol(SCZ_4)])

Control_1 = read.table('control-1-8-neuron.txt',header = T)
Control_2 = read.table('Control-9-16-neuron.txt',header = T)
Control_3 = read.table('Control-17-24-neuron.txt',header = T)
Control_4= read.table('Control-25-29-neuron.txt',header = T)
Control_count = data.frame(Control_1[,1],Control_1[,7:ncol(Control_1)],Control_2[,21:22],Control_2[,7:20],
                       Control_3[,7:ncol(Control_3)],Control_4[,7:ncol(Control_4)])
all_count = data.frame(SCZ_count,Control_count[,2:ncol(Control_count)])

colnames(all_count) = c('gene_id','SCZ_1_1','SCZ_1_2','SCZ_2_1','SCZ_2_2','SCZ_3_1','SCZ_3_2',
                        'SCZ_4_1','SCZ_4_2','SCZ_5_1','SCZ_5_2','SCZ_6_1','SCZ_6_2','SCZ_7_1','SCZ_7_2',
                        'SCZ_8_1','SCZ_8_2','SCZ_9_1','SCZ_9_2','SCZ_10_1','SCZ_10_2','SCZ_11_1','SCZ_11_2',
                        'SCZ_12_1','SCZ_12_2','SCZ_13_1','SCZ_13_2','SCZ_14_1','SCZ_14_2','SCZ_15_1','SCZ_15_2',
                        'SCZ_16_1','SCZ_16_2','SCZ_17_1','SCZ_17_2','SCZ_18_1','SCZ_18_2','SCZ_19_1','SCZ_19_2',
                        'SCZ_20_1','SCZ_20_2','SCZ_21_1','SCZ_21_2','SCZ_22_1','SCZ_22_2','SCZ_23_1','SCZ_23_2',
                        'SCZ_24_1','SCZ_24_2','SCZ_25_1','SCZ_25_2','SCZ_26_1','SCZ_26_2','SCZ_27_1','SCZ_27_2',
                        'SCZ_28_1','SCZ_28_2','SCZ_29_1','SCZ_29_2',
                        'Control_1_1','Control_1_2','Control_2_1','Control_2_2','Control_3_1','Control_3_2',
                        'Control_4_1','Control_4_2','Control_5_1','Control_5_2','Control_6_1','Control_6_2','Control_7_1','Control_7_2',
                        'Control_8_1','Control_8_2','Control_9_1','Control_9_2','Control_10_1','Control_10_2','Control_11_1','Control_11_2',
                        'Control_12_1','Control_12_2','Control_13_1','Control_13_2','Control_14_1','Control_14_2','Control_15_1','Control_15_2',
                        'Control_16_1','Control_16_2','Control_17_1','Control_17_2','Control_18_1','Control_18_2','Control_19_1','Control_19_2',
                        'Control_20_1','Control_20_2','Control_21_1','Control_21_2','Control_22_1','Control_22_2','Control_23_1','Control_23_2',
                        'Control_24_1','Control_24_2','Control_25_1','Control_25_2','Control_26_1','Control_26_2','Control_27_1','Control_27_2',
                        'Control_28_1','Control_28_2','Control_29_1','Control_29_2')

untrt_SCZ = data.frame(all_count[,1:31])
trt_SCZ = data.frame(all_count[,1],all_count[,32:59])
colnames(trt_SCZ)[1] = 'gene_id'

Control_set1 = data.frame(all_count[,1],all_count[,60:89])
colnames(Control_set1)[1] = 'gene_id'
Control_set2 = data.frame(all_count[,1],all_count[,90:117])
colnames(Control_set2)[1] = 'gene_id'

###############################SCZ vs. Control###############################################################
meta = factors
meta$condition_num = c(rep(1,58),rep(2,58))
SCZ_Control = all_count
# SCZ_Control = SCZ_Control[,-60]
row.names(SCZ_Control) = SCZ_Control[,1]
SCZ_Control = SCZ_Control[,-1]
miss_1 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,1:58]<20)) > 0.5*ncol(SCZ_Control[,1:58])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,59:116]<20)) > 0.5*ncol(SCZ_Control[,59:116])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
SCZ_Control = SCZ_Control[-miss,]
# t_SCZ_Control = t(SCZ_Control)
# umap_results <- umap(t_SCZ_Control,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
PCA_results = prcomp(SCZ_Control, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$total.reads,method = 'pearson')
cor_con_align = cor(meta$condition_num,meta$mapped.rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')
cor_con_dup = cor(meta$condition_num,meta$duplicate.rate,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$total.reads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$total.reads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$total.reads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$total.reads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$total.reads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$total.reads,method = 'pearson')

corr_uni_align_1 = cor(PCA$PC1,meta$mapped.rate,method = 'pearson')
corr_uni_align_2 = cor(PCA$PC2,meta$mapped.rate,method = 'pearson')
corr_uni_align_3 = cor(PCA$PC3,meta$mapped.rate,method = 'pearson')
corr_uni_align_4 = cor(PCA$PC4,meta$mapped.rate,method = 'pearson')
corr_uni_align_5 = cor(PCA$PC5,meta$mapped.rate,method = 'pearson')
corr_uni_align_6 = cor(PCA$PC6,meta$mapped.rate,method = 'pearson')

corr_uni_exon_1 = cor(PCA$PC1,meta$exon,method = 'pearson')
corr_uni_exon_2 = cor(PCA$PC2,meta$exon,method = 'pearson')
corr_uni_exon_3 = cor(PCA$PC3,meta$exon,method = 'pearson')
corr_uni_exon_4 = cor(PCA$PC4,meta$exon,method = 'pearson')
corr_uni_exon_5 = cor(PCA$PC5,meta$exon,method = 'pearson')
corr_uni_exon_6 = cor(PCA$PC6,meta$exon,method = 'pearson')

corr_uni_gender_1 = cor(PCA$PC1,meta$gender,method = 'pearson')
corr_uni_gender_2 = cor(PCA$PC2,meta$gender,method = 'pearson')
corr_uni_gender_3 = cor(PCA$PC3,meta$gender,method = 'pearson')
corr_uni_gender_4 = cor(PCA$PC4,meta$gender,method = 'pearson')
corr_uni_gender_5 = cor(PCA$PC5,meta$gender,method = 'pearson')
corr_uni_gender_6 = cor(PCA$PC6,meta$gender,method = 'pearson')

corr_uni_age_1 = cor(PCA$PC1,meta$age,method = 'pearson')
corr_uni_age_2 = cor(PCA$PC2,meta$age,method = 'pearson')
corr_uni_age_3 = cor(PCA$PC3,meta$age,method = 'pearson')
corr_uni_age_4 = cor(PCA$PC4,meta$age,method = 'pearson')
corr_uni_age_5 = cor(PCA$PC5,meta$age,method = 'pearson')
corr_uni_age_6 = cor(PCA$PC6,meta$age,method = 'pearson')

corr_uni_PMI_1 = cor(PCA$PC1,meta$PMI,method = 'pearson')
corr_uni_PMI_2 = cor(PCA$PC2,meta$PMI,method = 'pearson')
corr_uni_PMI_3 = cor(PCA$PC3,meta$PMI,method = 'pearson')
corr_uni_PMI_4 = cor(PCA$PC4,meta$PMI,method = 'pearson')
corr_uni_PMI_5 = cor(PCA$PC5,meta$PMI,method = 'pearson')
corr_uni_PMI_6 = cor(PCA$PC6,meta$PMI,method = 'pearson')

corr_uni_dup_1 = cor(PCA$PC1,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_2 = cor(PCA$PC2,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_3 = cor(PCA$PC3,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_4 = cor(PCA$PC4,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_5 = cor(PCA$PC5,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_6 = cor(PCA$PC6,meta$duplicate.rate,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$total.reads
group_list3 = meta$mapped.rate
group_list4 = meta$exon
group_list5 = meta$gender
group_list6 = meta$age
group_list7 = meta$PMI
group_list8 = meta$duplicate.rate

colData=data.frame(row.names = colnames(SCZ_Control),
                   condition=group_list1,
                   total=group_list2,
                   align=group_list3,
                   exon=group_list4,
                   gender=group_list5,
                   age = group_list6,
                   PMI = group_list7,
                   dup = group_list8)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$align=as.numeric(colData$align)
colData$exon=as.numeric(colData$exon)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=as.numeric(colData$PMI)
colData$dup=as.numeric(colData$dup)

dds=DESeqDataSetFromMatrix(countData = SCZ_Control,
                           colData = colData,
                           design = ~ condition + age + align + dup + exon)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=SCZ_Control[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/neuron/")
all_glia = SCZ_Control[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_SCZ_Control_neuron.csv')
write.csv(DEG_matrix,'Differential_SCZ_Control_neuron.csv')



###############################untrt_SCZ vs. Control_1###############################################################
meta = untrt_Control1_factors
meta$condition_num = c(rep(1,30),rep(2,30))
SCZ_Control = data.frame(untrt_SCZ,Control_set1)
SCZ_Control = SCZ_Control[,-32]
row.names(SCZ_Control) = SCZ_Control[,1]
SCZ_Control = SCZ_Control[,-1]
miss_1 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,1:30]<20)) > 0.5*ncol(SCZ_Control[,1:30])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,31:60]<20)) > 0.5*ncol(SCZ_Control[,31:60])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
SCZ_Control = SCZ_Control[-miss,]
# t_SCZ_Control = t(SCZ_Control)
# umap_results <- umap(t_SCZ_Control,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
PCA_results = prcomp(SCZ_Control, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$total.reads,method = 'pearson')
cor_con_align = cor(meta$condition_num,meta$mapped.rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')
cor_con_dup = cor(meta$condition_num,meta$duplicate.rate,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$total.reads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$total.reads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$total.reads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$total.reads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$total.reads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$total.reads,method = 'pearson')

corr_uni_align_1 = cor(PCA$PC1,meta$mapped.rate,method = 'pearson')
corr_uni_align_2 = cor(PCA$PC2,meta$mapped.rate,method = 'pearson')
corr_uni_align_3 = cor(PCA$PC3,meta$mapped.rate,method = 'pearson')
corr_uni_align_4 = cor(PCA$PC4,meta$mapped.rate,method = 'pearson')
corr_uni_align_5 = cor(PCA$PC5,meta$mapped.rate,method = 'pearson')
corr_uni_align_6 = cor(PCA$PC6,meta$mapped.rate,method = 'pearson')

corr_uni_exon_1 = cor(PCA$PC1,meta$exon,method = 'pearson')
corr_uni_exon_2 = cor(PCA$PC2,meta$exon,method = 'pearson')
corr_uni_exon_3 = cor(PCA$PC3,meta$exon,method = 'pearson')
corr_uni_exon_4 = cor(PCA$PC4,meta$exon,method = 'pearson')
corr_uni_exon_5 = cor(PCA$PC5,meta$exon,method = 'pearson')
corr_uni_exon_6 = cor(PCA$PC6,meta$exon,method = 'pearson')

corr_uni_gender_1 = cor(PCA$PC1,meta$gender,method = 'pearson')
corr_uni_gender_2 = cor(PCA$PC2,meta$gender,method = 'pearson')
corr_uni_gender_3 = cor(PCA$PC3,meta$gender,method = 'pearson')
corr_uni_gender_4 = cor(PCA$PC4,meta$gender,method = 'pearson')
corr_uni_gender_5 = cor(PCA$PC5,meta$gender,method = 'pearson')
corr_uni_gender_6 = cor(PCA$PC6,meta$gender,method = 'pearson')

corr_uni_age_1 = cor(PCA$PC1,meta$age,method = 'pearson')
corr_uni_age_2 = cor(PCA$PC2,meta$age,method = 'pearson')
corr_uni_age_3 = cor(PCA$PC3,meta$age,method = 'pearson')
corr_uni_age_4 = cor(PCA$PC4,meta$age,method = 'pearson')
corr_uni_age_5 = cor(PCA$PC5,meta$age,method = 'pearson')
corr_uni_age_6 = cor(PCA$PC6,meta$age,method = 'pearson')

corr_uni_PMI_1 = cor(PCA$PC1,meta$PMI,method = 'pearson')
corr_uni_PMI_2 = cor(PCA$PC2,meta$PMI,method = 'pearson')
corr_uni_PMI_3 = cor(PCA$PC3,meta$PMI,method = 'pearson')
corr_uni_PMI_4 = cor(PCA$PC4,meta$PMI,method = 'pearson')
corr_uni_PMI_5 = cor(PCA$PC5,meta$PMI,method = 'pearson')
corr_uni_PMI_6 = cor(PCA$PC6,meta$PMI,method = 'pearson')

corr_uni_dup_1 = cor(PCA$PC1,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_2 = cor(PCA$PC2,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_3 = cor(PCA$PC3,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_4 = cor(PCA$PC4,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_5 = cor(PCA$PC5,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_6 = cor(PCA$PC6,meta$duplicate.rate,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$total.reads
group_list3 = meta$mapped.rate
group_list4 = meta$exon
group_list5 = meta$gender
group_list6 = meta$age
group_list7 = meta$PMI
group_list8 = meta$duplicate.rate

colData=data.frame(row.names = colnames(SCZ_Control),
                   condition=group_list1,
                   total=group_list2,
                   align=group_list3,
                   exon=group_list4,
                   gender=group_list5,
                   age = group_list6,
                   PMI = group_list7,
                   dup = group_list8)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$align=as.numeric(colData$align)
colData$exon=as.numeric(colData$exon)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=as.numeric(colData$PMI)
colData$dup=as.numeric(colData$dup)

dds=DESeqDataSetFromMatrix(countData = SCZ_Control,
                           colData = colData,
                           design = ~ condition + age + align + dup + exon + PMI)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=SCZ_Control[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/neuron/")
all_glia = SCZ_Control[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_untrt_SCZ_Control_1_neuron.csv')
write.csv(DEG_matrix,'Differential_untrt_SCZ_Control_1_neuron.csv')


###############################trt_SCZ vs. Control_2###############################################################
meta = trt_Control2_factors
meta$condition_num = c(rep(1,28),rep(2,28))
SCZ_Control = data.frame(trt_SCZ,Control_set2)
SCZ_Control = SCZ_Control[,-30]
row.names(SCZ_Control) = SCZ_Control[,1]
SCZ_Control = SCZ_Control[,-1]
miss_1 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,1:28]<20)) > 0.5*ncol(SCZ_Control[,1:28])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,29:56]<20)) > 0.5*ncol(SCZ_Control[,29:56])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
SCZ_Control = SCZ_Control[-miss,]
# t_SCZ_Control = t(SCZ_Control)
# umap_results <- umap(t_SCZ_Control,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
PCA_results = prcomp(SCZ_Control, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$total.reads,method = 'pearson')
cor_con_align = cor(meta$condition_num,meta$mapped.rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')
cor_con_dup = cor(meta$condition_num,meta$duplicate.rate,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$total.reads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$total.reads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$total.reads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$total.reads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$total.reads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$total.reads,method = 'pearson')

corr_uni_align_1 = cor(PCA$PC1,meta$mapped.rate,method = 'pearson')
corr_uni_align_2 = cor(PCA$PC2,meta$mapped.rate,method = 'pearson')
corr_uni_align_3 = cor(PCA$PC3,meta$mapped.rate,method = 'pearson')
corr_uni_align_4 = cor(PCA$PC4,meta$mapped.rate,method = 'pearson')
corr_uni_align_5 = cor(PCA$PC5,meta$mapped.rate,method = 'pearson')
corr_uni_align_6 = cor(PCA$PC6,meta$mapped.rate,method = 'pearson')

corr_uni_exon_1 = cor(PCA$PC1,meta$exon,method = 'pearson')
corr_uni_exon_2 = cor(PCA$PC2,meta$exon,method = 'pearson')
corr_uni_exon_3 = cor(PCA$PC3,meta$exon,method = 'pearson')
corr_uni_exon_4 = cor(PCA$PC4,meta$exon,method = 'pearson')
corr_uni_exon_5 = cor(PCA$PC5,meta$exon,method = 'pearson')
corr_uni_exon_6 = cor(PCA$PC6,meta$exon,method = 'pearson')

corr_uni_gender_1 = cor(PCA$PC1,meta$gender,method = 'pearson')
corr_uni_gender_2 = cor(PCA$PC2,meta$gender,method = 'pearson')
corr_uni_gender_3 = cor(PCA$PC3,meta$gender,method = 'pearson')
corr_uni_gender_4 = cor(PCA$PC4,meta$gender,method = 'pearson')
corr_uni_gender_5 = cor(PCA$PC5,meta$gender,method = 'pearson')
corr_uni_gender_6 = cor(PCA$PC6,meta$gender,method = 'pearson')

corr_uni_age_1 = cor(PCA$PC1,meta$age,method = 'pearson')
corr_uni_age_2 = cor(PCA$PC2,meta$age,method = 'pearson')
corr_uni_age_3 = cor(PCA$PC3,meta$age,method = 'pearson')
corr_uni_age_4 = cor(PCA$PC4,meta$age,method = 'pearson')
corr_uni_age_5 = cor(PCA$PC5,meta$age,method = 'pearson')
corr_uni_age_6 = cor(PCA$PC6,meta$age,method = 'pearson')

corr_uni_PMI_1 = cor(PCA$PC1,meta$PMI,method = 'pearson')
corr_uni_PMI_2 = cor(PCA$PC2,meta$PMI,method = 'pearson')
corr_uni_PMI_3 = cor(PCA$PC3,meta$PMI,method = 'pearson')
corr_uni_PMI_4 = cor(PCA$PC4,meta$PMI,method = 'pearson')
corr_uni_PMI_5 = cor(PCA$PC5,meta$PMI,method = 'pearson')
corr_uni_PMI_6 = cor(PCA$PC6,meta$PMI,method = 'pearson')

corr_uni_dup_1 = cor(PCA$PC1,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_2 = cor(PCA$PC2,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_3 = cor(PCA$PC3,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_4 = cor(PCA$PC4,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_5 = cor(PCA$PC5,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_6 = cor(PCA$PC6,meta$duplicate.rate,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$total.reads
group_list3 = meta$mapped.rate
group_list4 = meta$exon
group_list5 = meta$gender
group_list6 = meta$age
group_list7 = meta$PMI
group_list8 = meta$duplicate.rate

colData=data.frame(row.names = colnames(SCZ_Control),
                   condition=group_list1,
                   total=group_list2,
                   align=group_list3,
                   exon=group_list4,
                   gender=group_list5,
                   age = group_list6,
                   PMI = group_list7,
                   dup = group_list8)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$align=as.numeric(colData$align)
colData$exon=as.numeric(colData$exon)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=as.numeric(colData$PMI)
colData$dup=as.numeric(colData$dup)

dds=DESeqDataSetFromMatrix(countData = SCZ_Control,
                           colData = colData,
                           design = ~ condition + age + align + dup + exon + total)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=SCZ_Control[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/neuron/")
all_glia = SCZ_Control[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_trt_SCZ_Control_2_neuron.csv')
write.csv(DEG_matrix,'Differential_trt_SCZ_Control_2_neuron.csv')

###################glia#############################
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/")
factors = read.csv('metadata_RNA_glia.csv')
row.names(factors) = factors[,1]
factors = factors[,-1]
factors$total.reads = round((factors$total.reads/1000000),2)
untrt_Control1_factors = bind_rows(factors[1:30,],factors[59:88,])
trt_Control2_factors = bind_rows(factors[31:58,],factors[89:116,])

setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/glia/")
SCZ_1 = read.table('SCZ-1-8-glia.txt',header = T)
SCZ_2 = read.table('SCZ-9-16-glia.txt',header = T)
SCZ_3 = read.table('SCZ-17-24-glia.txt',header = T)
SCZ_4= read.table('SCZ-25-29-glia.txt',header = T)
SCZ_count = data.frame(SCZ_1[,1],SCZ_1[,7:ncol(SCZ_1)],SCZ_2[,21:22],SCZ_2[,7:20],
                       SCZ_3[,7:ncol(SCZ_3)],SCZ_4[,7:ncol(SCZ_4)])

Control_1 = read.table('control-1-8-glia.txt',header = T)
Control_2 = read.table('control-9-16-glia.txt',header = T)
Control_3 = read.table('control-17-24-glia.txt',header = T)
Control_4= read.table('control-25-29-glia.txt',header = T)
Control_count = data.frame(Control_1[,1],Control_1[,7:ncol(Control_1)],Control_2[,21:22],Control_2[,7:20],
                           Control_3[,7:ncol(Control_3)],Control_4[,7:ncol(Control_4)])
all_count = data.frame(SCZ_count,Control_count[,2:ncol(Control_count)])

colnames(all_count) = c('gene_id','SCZ_1_1','SCZ_1_2','SCZ_2_1','SCZ_2_2','SCZ_3_1','SCZ_3_2',
                        'SCZ_4_1','SCZ_4_2','SCZ_5_1','SCZ_5_2','SCZ_6_1','SCZ_6_2','SCZ_7_1','SCZ_7_2',
                        'SCZ_8_1','SCZ_8_2','SCZ_9_1','SCZ_9_2','SCZ_10_1','SCZ_10_2','SCZ_11_1','SCZ_11_2',
                        'SCZ_12_1','SCZ_12_2','SCZ_13_1','SCZ_13_2','SCZ_14_1','SCZ_14_2','SCZ_15_1','SCZ_15_2',
                        'SCZ_16_1','SCZ_16_2','SCZ_17_1','SCZ_17_2','SCZ_18_1','SCZ_18_2','SCZ_19_1','SCZ_19_2',
                        'SCZ_20_1','SCZ_20_2','SCZ_21_1','SCZ_21_2','SCZ_22_1','SCZ_22_2','SCZ_23_1','SCZ_23_2',
                        'SCZ_24_1','SCZ_24_2','SCZ_25_1','SCZ_25_2','SCZ_26_1','SCZ_26_2','SCZ_27_1','SCZ_27_2',
                        'SCZ_28_1','SCZ_28_2','SCZ_29_1','SCZ_29_2',
                        'Control_1_1','Control_1_2','Control_2_1','Control_2_2','Control_3_1','Control_3_2',
                        'Control_4_1','Control_4_2','Control_5_1','Control_5_2','Control_6_1','Control_6_2','Control_7_1','Control_7_2',
                        'Control_8_1','Control_8_2','Control_9_1','Control_9_2','Control_10_1','Control_10_2','Control_11_1','Control_11_2',
                        'Control_12_1','Control_12_2','Control_13_1','Control_13_2','Control_14_1','Control_14_2','Control_15_1','Control_15_2',
                        'Control_16_1','Control_16_2','Control_17_1','Control_17_2','Control_18_1','Control_18_2','Control_19_1','Control_19_2',
                        'Control_20_1','Control_20_2','Control_21_1','Control_21_2','Control_22_1','Control_22_2','Control_23_1','Control_23_2',
                        'Control_24_1','Control_24_2','Control_25_1','Control_25_2','Control_26_1','Control_26_2','Control_27_1','Control_27_2',
                        'Control_28_1','Control_28_2','Control_29_1','Control_29_2')

untrt_SCZ = data.frame(all_count[,1:31])
trt_SCZ = data.frame(all_count[,1],all_count[,32:59])
colnames(trt_SCZ)[1] = 'gene_id'

Control_set1 = data.frame(all_count[,1],all_count[,60:89])
colnames(Control_set1)[1] = 'gene_id'
Control_set2 = data.frame(all_count[,1],all_count[,90:117])
colnames(Control_set2)[1] = 'gene_id'

###############################SCZ vs. Control###############################################################
meta = factors
meta$condition_num = c(rep(1,58),rep(2,58))
SCZ_Control = all_count
# SCZ_Control = SCZ_Control[,-60]
row.names(SCZ_Control) = SCZ_Control[,1]
SCZ_Control = SCZ_Control[,-1]
miss_1 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,1:58]<20)) > 0.5*ncol(SCZ_Control[,1:58])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,59:116]<20)) > 0.5*ncol(SCZ_Control[,59:116])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
SCZ_Control = SCZ_Control[-miss,]
# t_SCZ_Control = t(SCZ_Control)
# umap_results <- umap(t_SCZ_Control,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
PCA_results = prcomp(SCZ_Control, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$total.reads,method = 'pearson')
cor_con_align = cor(meta$condition_num,meta$mapped.rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')
cor_con_dup = cor(meta$condition_num,meta$duplicate.rate,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$total.reads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$total.reads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$total.reads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$total.reads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$total.reads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$total.reads,method = 'pearson')

corr_uni_align_1 = cor(PCA$PC1,meta$mapped.rate,method = 'pearson')
corr_uni_align_2 = cor(PCA$PC2,meta$mapped.rate,method = 'pearson')
corr_uni_align_3 = cor(PCA$PC3,meta$mapped.rate,method = 'pearson')
corr_uni_align_4 = cor(PCA$PC4,meta$mapped.rate,method = 'pearson')
corr_uni_align_5 = cor(PCA$PC5,meta$mapped.rate,method = 'pearson')
corr_uni_align_6 = cor(PCA$PC6,meta$mapped.rate,method = 'pearson')

corr_uni_exon_1 = cor(PCA$PC1,meta$exon,method = 'pearson')
corr_uni_exon_2 = cor(PCA$PC2,meta$exon,method = 'pearson')
corr_uni_exon_3 = cor(PCA$PC3,meta$exon,method = 'pearson')
corr_uni_exon_4 = cor(PCA$PC4,meta$exon,method = 'pearson')
corr_uni_exon_5 = cor(PCA$PC5,meta$exon,method = 'pearson')
corr_uni_exon_6 = cor(PCA$PC6,meta$exon,method = 'pearson')

corr_uni_gender_1 = cor(PCA$PC1,meta$gender,method = 'pearson')
corr_uni_gender_2 = cor(PCA$PC2,meta$gender,method = 'pearson')
corr_uni_gender_3 = cor(PCA$PC3,meta$gender,method = 'pearson')
corr_uni_gender_4 = cor(PCA$PC4,meta$gender,method = 'pearson')
corr_uni_gender_5 = cor(PCA$PC5,meta$gender,method = 'pearson')
corr_uni_gender_6 = cor(PCA$PC6,meta$gender,method = 'pearson')

corr_uni_age_1 = cor(PCA$PC1,meta$age,method = 'pearson')
corr_uni_age_2 = cor(PCA$PC2,meta$age,method = 'pearson')
corr_uni_age_3 = cor(PCA$PC3,meta$age,method = 'pearson')
corr_uni_age_4 = cor(PCA$PC4,meta$age,method = 'pearson')
corr_uni_age_5 = cor(PCA$PC5,meta$age,method = 'pearson')
corr_uni_age_6 = cor(PCA$PC6,meta$age,method = 'pearson')

corr_uni_PMI_1 = cor(PCA$PC1,meta$PMI,method = 'pearson')
corr_uni_PMI_2 = cor(PCA$PC2,meta$PMI,method = 'pearson')
corr_uni_PMI_3 = cor(PCA$PC3,meta$PMI,method = 'pearson')
corr_uni_PMI_4 = cor(PCA$PC4,meta$PMI,method = 'pearson')
corr_uni_PMI_5 = cor(PCA$PC5,meta$PMI,method = 'pearson')
corr_uni_PMI_6 = cor(PCA$PC6,meta$PMI,method = 'pearson')

corr_uni_dup_1 = cor(PCA$PC1,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_2 = cor(PCA$PC2,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_3 = cor(PCA$PC3,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_4 = cor(PCA$PC4,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_5 = cor(PCA$PC5,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_6 = cor(PCA$PC6,meta$duplicate.rate,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$total.reads
group_list3 = meta$mapped.rate
group_list4 = meta$exon
group_list5 = meta$gender
group_list6 = meta$age
group_list7 = meta$PMI
group_list8 = meta$duplicate.rate

colData=data.frame(row.names = colnames(SCZ_Control),
                   condition=group_list1,
                   total=group_list2,
                   align=group_list3,
                   exon=group_list4,
                   gender=group_list5,
                   age = group_list6,
                   PMI = group_list7,
                   dup = group_list8)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$align=as.numeric(colData$align)
colData$exon=as.numeric(colData$exon)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=as.numeric(colData$PMI)
colData$dup=as.numeric(colData$dup)

dds=DESeqDataSetFromMatrix(countData = SCZ_Control,
                           colData = colData,
                           design = ~ condition)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=SCZ_Control[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/glia/")
all_glia = SCZ_Control[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_SCZ_Control_glia.csv')
write.csv(DEG_matrix,'Differential_SCZ_Control_glia.csv')


###############################untrt_SCZ vs. Control_1###############################################################
meta = untrt_Control1_factors
meta$condition_num = c(rep(1,30),rep(2,30))
SCZ_Control = data.frame(untrt_SCZ,Control_set1)
SCZ_Control = SCZ_Control[,-32]
row.names(SCZ_Control) = SCZ_Control[,1]
SCZ_Control = SCZ_Control[,-1]
miss_1 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,1:30]<20)) > 0.5*ncol(SCZ_Control[,1:30])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,31:60]<20)) > 0.5*ncol(SCZ_Control[,31:60])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
SCZ_Control = SCZ_Control[-miss,]
# t_SCZ_Control = t(SCZ_Control)
# umap_results <- umap(t_SCZ_Control,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
PCA_results = prcomp(SCZ_Control, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$total.reads,method = 'pearson')
cor_con_align = cor(meta$condition_num,meta$mapped.rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')
cor_con_dup = cor(meta$condition_num,meta$duplicate.rate,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$total.reads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$total.reads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$total.reads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$total.reads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$total.reads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$total.reads,method = 'pearson')

corr_uni_align_1 = cor(PCA$PC1,meta$mapped.rate,method = 'pearson')
corr_uni_align_2 = cor(PCA$PC2,meta$mapped.rate,method = 'pearson')
corr_uni_align_3 = cor(PCA$PC3,meta$mapped.rate,method = 'pearson')
corr_uni_align_4 = cor(PCA$PC4,meta$mapped.rate,method = 'pearson')
corr_uni_align_5 = cor(PCA$PC5,meta$mapped.rate,method = 'pearson')
corr_uni_align_6 = cor(PCA$PC6,meta$mapped.rate,method = 'pearson')

corr_uni_exon_1 = cor(PCA$PC1,meta$exon,method = 'pearson')
corr_uni_exon_2 = cor(PCA$PC2,meta$exon,method = 'pearson')
corr_uni_exon_3 = cor(PCA$PC3,meta$exon,method = 'pearson')
corr_uni_exon_4 = cor(PCA$PC4,meta$exon,method = 'pearson')
corr_uni_exon_5 = cor(PCA$PC5,meta$exon,method = 'pearson')
corr_uni_exon_6 = cor(PCA$PC6,meta$exon,method = 'pearson')

corr_uni_gender_1 = cor(PCA$PC1,meta$gender,method = 'pearson')
corr_uni_gender_2 = cor(PCA$PC2,meta$gender,method = 'pearson')
corr_uni_gender_3 = cor(PCA$PC3,meta$gender,method = 'pearson')
corr_uni_gender_4 = cor(PCA$PC4,meta$gender,method = 'pearson')
corr_uni_gender_5 = cor(PCA$PC5,meta$gender,method = 'pearson')
corr_uni_gender_6 = cor(PCA$PC6,meta$gender,method = 'pearson')

corr_uni_age_1 = cor(PCA$PC1,meta$age,method = 'pearson')
corr_uni_age_2 = cor(PCA$PC2,meta$age,method = 'pearson')
corr_uni_age_3 = cor(PCA$PC3,meta$age,method = 'pearson')
corr_uni_age_4 = cor(PCA$PC4,meta$age,method = 'pearson')
corr_uni_age_5 = cor(PCA$PC5,meta$age,method = 'pearson')
corr_uni_age_6 = cor(PCA$PC6,meta$age,method = 'pearson')

corr_uni_PMI_1 = cor(PCA$PC1,meta$PMI,method = 'pearson')
corr_uni_PMI_2 = cor(PCA$PC2,meta$PMI,method = 'pearson')
corr_uni_PMI_3 = cor(PCA$PC3,meta$PMI,method = 'pearson')
corr_uni_PMI_4 = cor(PCA$PC4,meta$PMI,method = 'pearson')
corr_uni_PMI_5 = cor(PCA$PC5,meta$PMI,method = 'pearson')
corr_uni_PMI_6 = cor(PCA$PC6,meta$PMI,method = 'pearson')

corr_uni_dup_1 = cor(PCA$PC1,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_2 = cor(PCA$PC2,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_3 = cor(PCA$PC3,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_4 = cor(PCA$PC4,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_5 = cor(PCA$PC5,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_6 = cor(PCA$PC6,meta$duplicate.rate,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$total.reads
group_list3 = meta$mapped.rate
group_list4 = meta$exon
group_list5 = meta$gender
group_list6 = meta$age
group_list7 = meta$PMI
group_list8 = meta$duplicate.rate

colData=data.frame(row.names = colnames(SCZ_Control),
                   condition=group_list1,
                   total=group_list2,
                   align=group_list3,
                   exon=group_list4,
                   gender=group_list5,
                   age = group_list6,
                   PMI = group_list7,
                   dup = group_list8)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$align=as.numeric(colData$align)
colData$exon=as.numeric(colData$exon)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=as.numeric(colData$PMI)
colData$dup=as.numeric(colData$dup)

dds=DESeqDataSetFromMatrix(countData = SCZ_Control,
                           colData = colData,
                           design = ~ condition + exon)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=SCZ_Control[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/glia/")
all_glia = SCZ_Control[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_untrt_SCZ_Control_1_glia.csv')
write.csv(DEG_matrix,'Differential_untrt_SCZ_Control_1_glia.csv')


###############################trt_SCZ vs. Control_2###############################################################
meta = trt_Control2_factors
meta$condition_num = c(rep(1,28),rep(2,28))
SCZ_Control = data.frame(trt_SCZ,Control_set2)
SCZ_Control = SCZ_Control[,-30]
row.names(SCZ_Control) = SCZ_Control[,1]
SCZ_Control = SCZ_Control[,-1]
miss_1 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,1:28]<20)) > 0.5*ncol(SCZ_Control[,1:28])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(SCZ_Control)) {
  if(length(which(SCZ_Control[i,29:56]<20)) > 0.5*ncol(SCZ_Control[,29:56])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
SCZ_Control = SCZ_Control[-miss,]
# t_SCZ_Control = t(SCZ_Control)
# umap_results <- umap(t_SCZ_Control,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
PCA_results = prcomp(SCZ_Control, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$total.reads,method = 'pearson')
cor_con_align = cor(meta$condition_num,meta$mapped.rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')
cor_con_dup = cor(meta$condition_num,meta$duplicate.rate,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$total.reads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$total.reads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$total.reads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$total.reads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$total.reads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$total.reads,method = 'pearson')

corr_uni_align_1 = cor(PCA$PC1,meta$mapped.rate,method = 'pearson')
corr_uni_align_2 = cor(PCA$PC2,meta$mapped.rate,method = 'pearson')
corr_uni_align_3 = cor(PCA$PC3,meta$mapped.rate,method = 'pearson')
corr_uni_align_4 = cor(PCA$PC4,meta$mapped.rate,method = 'pearson')
corr_uni_align_5 = cor(PCA$PC5,meta$mapped.rate,method = 'pearson')
corr_uni_align_6 = cor(PCA$PC6,meta$mapped.rate,method = 'pearson')

corr_uni_exon_1 = cor(PCA$PC1,meta$exon,method = 'pearson')
corr_uni_exon_2 = cor(PCA$PC2,meta$exon,method = 'pearson')
corr_uni_exon_3 = cor(PCA$PC3,meta$exon,method = 'pearson')
corr_uni_exon_4 = cor(PCA$PC4,meta$exon,method = 'pearson')
corr_uni_exon_5 = cor(PCA$PC5,meta$exon,method = 'pearson')
corr_uni_exon_6 = cor(PCA$PC6,meta$exon,method = 'pearson')

corr_uni_gender_1 = cor(PCA$PC1,meta$gender,method = 'pearson')
corr_uni_gender_2 = cor(PCA$PC2,meta$gender,method = 'pearson')
corr_uni_gender_3 = cor(PCA$PC3,meta$gender,method = 'pearson')
corr_uni_gender_4 = cor(PCA$PC4,meta$gender,method = 'pearson')
corr_uni_gender_5 = cor(PCA$PC5,meta$gender,method = 'pearson')
corr_uni_gender_6 = cor(PCA$PC6,meta$gender,method = 'pearson')

corr_uni_age_1 = cor(PCA$PC1,meta$age,method = 'pearson')
corr_uni_age_2 = cor(PCA$PC2,meta$age,method = 'pearson')
corr_uni_age_3 = cor(PCA$PC3,meta$age,method = 'pearson')
corr_uni_age_4 = cor(PCA$PC4,meta$age,method = 'pearson')
corr_uni_age_5 = cor(PCA$PC5,meta$age,method = 'pearson')
corr_uni_age_6 = cor(PCA$PC6,meta$age,method = 'pearson')

corr_uni_PMI_1 = cor(PCA$PC1,meta$PMI,method = 'pearson')
corr_uni_PMI_2 = cor(PCA$PC2,meta$PMI,method = 'pearson')
corr_uni_PMI_3 = cor(PCA$PC3,meta$PMI,method = 'pearson')
corr_uni_PMI_4 = cor(PCA$PC4,meta$PMI,method = 'pearson')
corr_uni_PMI_5 = cor(PCA$PC5,meta$PMI,method = 'pearson')
corr_uni_PMI_6 = cor(PCA$PC6,meta$PMI,method = 'pearson')

corr_uni_dup_1 = cor(PCA$PC1,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_2 = cor(PCA$PC2,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_3 = cor(PCA$PC3,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_4 = cor(PCA$PC4,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_5 = cor(PCA$PC5,meta$duplicate.rate,method = 'pearson')
corr_uni_dup_6 = cor(PCA$PC6,meta$duplicate.rate,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$total.reads
group_list3 = meta$mapped.rate
group_list4 = meta$exon
group_list5 = meta$gender
group_list6 = meta$age
group_list7 = meta$PMI
group_list8 = meta$duplicate.rate

colData=data.frame(row.names = colnames(SCZ_Control),
                   condition=group_list1,
                   total=group_list2,
                   align=group_list3,
                   exon=group_list4,
                   gender=group_list5,
                   age = group_list6,
                   PMI = group_list7,
                   dup = group_list8)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$align=as.numeric(colData$align)
colData$exon=as.numeric(colData$exon)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=as.numeric(colData$PMI)
colData$dup=as.numeric(colData$dup)

dds=DESeqDataSetFromMatrix(countData = SCZ_Control,
                           colData = colData,
                           design = ~ condition + align + total)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=SCZ_Control[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/glia/")
all_glia = SCZ_Control[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_trt_SCZ_Control_2_glia.csv')
write.csv(DEG_matrix,'Differential_trt_SCZ_Control_2_glia.csv')

