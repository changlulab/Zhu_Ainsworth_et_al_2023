##############################neuron###########
############preparation##############
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/cofactors/")
factors = read.csv('enhancer_SCZ_Control_neuron.csv')
row.names(factors) = factors[,1]
factors = factors[,-1]
factors$totalreads = round((factors$totalreads/1000000),2)
factors$uniquereads = round((factors$uniquereads/1000000),2)
factors$peaks = round((factors$peaks/1000),2)
factors$NSC = round(factors$NSC,2)
factors$RSC = round(factors$RSC,2)
factors$PBC = round(factors$PBC,2)
factors_all = factors
factors_untrt = bind_rows(factors[1:30,],factors[59:88,])
factors_trt = bind_rows(factors[31:58,],factors[89:116,])
###############################SCZ vs. Control###############################################################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/count files/")
raw_def = read.csv('SCZ_Control_neuron_enhancer_count.csv')
raw_def = all_enhancer
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
meta = factors_all
meta$condition = c(rep('SCZ',58),rep('Control',58))
meta$condition_num = c(rep(1,58),rep(2,58))
group_all = def
miss_1 <- c()
for (i in 1:nrow(group_all)) {
  if(length(which(group_all[i,1:58]<20)) > 0.5*ncol(group_all[,1:58])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_all)) {
  if(length(which(group_all[i,59:116]<20)) > 0.5*ncol(group_all[,59:116])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_all = group_all[-miss,]
# t_group_all = t(group_all)
# umap_results <- umap(t_group_all,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
PCA_results = prcomp(group_all, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$totalreads,method = 'pearson')
cor_con_unique = cor(meta$condition_num,meta$uniquereads,method = 'pearson')
cor_con_peak = cor(meta$condition_num,meta$peaks,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_PBC = cor(meta$condition_num,meta$PBC,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$totalreads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$totalreads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$totalreads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$totalreads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$totalreads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$totalreads,method = 'pearson')

corr_uni_unique_1 = cor(PCA$PC1,meta$uniquereads,method = 'pearson')
corr_uni_unique_2 = cor(PCA$PC2,meta$uniquereads,method = 'pearson')
corr_uni_unique_3 = cor(PCA$PC3,meta$uniquereads,method = 'pearson')
corr_uni_unique_4 = cor(PCA$PC4,meta$uniquereads,method = 'pearson')
corr_uni_unique_5 = cor(PCA$PC5,meta$uniquereads,method = 'pearson')
corr_uni_unique_6 = cor(PCA$PC6,meta$uniquereads,method = 'pearson')

corr_uni_peaks_1 = cor(PCA$PC1,meta$peaks,method = 'pearson')
corr_uni_peaks_2 = cor(PCA$PC2,meta$peaks,method = 'pearson')
corr_uni_peaks_3 = cor(PCA$PC3,meta$peaks,method = 'pearson')
corr_uni_peaks_4 = cor(PCA$PC4,meta$peaks,method = 'pearson')
corr_uni_peaks_5 = cor(PCA$PC5,meta$peaks,method = 'pearson')
corr_uni_peaks_6 = cor(PCA$PC6,meta$peaks,method = 'pearson')

corr_uni_FRiP_1 = cor(PCA$PC1,meta$FRiP,method = 'pearson')
corr_uni_FRiP_2 = cor(PCA$PC2,meta$FRiP,method = 'pearson')
corr_uni_FRiP_3 = cor(PCA$PC3,meta$FRiP,method = 'pearson')
corr_uni_FRiP_4 = cor(PCA$PC4,meta$FRiP,method = 'pearson')
corr_uni_FRiP_5 = cor(PCA$PC5,meta$FRiP,method = 'pearson')
corr_uni_FRiP_6 = cor(PCA$PC6,meta$FRiP,method = 'pearson')

corr_uni_NSC_1 = cor(PCA$PC1,meta$NSC,method = 'pearson')
corr_uni_NSC_2 = cor(PCA$PC2,meta$NSC,method = 'pearson')
corr_uni_NSC_3 = cor(PCA$PC3,meta$NSC,method = 'pearson')
corr_uni_NSC_4 = cor(PCA$PC4,meta$NSC,method = 'pearson')
corr_uni_NSC_5 = cor(PCA$PC5,meta$NSC,method = 'pearson')
corr_uni_NSC_6 = cor(PCA$PC6,meta$NSC,method = 'pearson')

corr_uni_RSC_1 = cor(PCA$PC1,meta$RSC,method = 'pearson')
corr_uni_RSC_2 = cor(PCA$PC2,meta$RSC,method = 'pearson')
corr_uni_RSC_3 = cor(PCA$PC3,meta$RSC,method = 'pearson')
corr_uni_RSC_4 = cor(PCA$PC4,meta$RSC,method = 'pearson')
corr_uni_RSC_5 = cor(PCA$PC5,meta$RSC,method = 'pearson')
corr_uni_RSC_6 = cor(PCA$PC6,meta$RSC,method = 'pearson')

corr_uni_PBC_1 = cor(PCA$PC1,meta$PBC,method = 'pearson')
corr_uni_PBC_2 = cor(PCA$PC2,meta$PBC,method = 'pearson')
corr_uni_PBC_3 = cor(PCA$PC3,meta$PBC,method = 'pearson')
corr_uni_PBC_4 = cor(PCA$PC4,meta$PBC,method = 'pearson')
corr_uni_PBC_5 = cor(PCA$PC5,meta$PBC,method = 'pearson')
corr_uni_PBC_6 = cor(PCA$PC6,meta$PBC,method = 'pearson')

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

group_list1 = meta$condition
group_list2 = meta$totalreads
group_list3 = meta$uniquereads
group_list4 = meta$peaks
group_list5 = meta$FRiP
group_list6 = meta$NSC
group_list7 = meta$RSC
group_list8 = meta$PBC
group_list9 = meta$gender
group_list10 = meta$age
group_list11 = meta$PMI

colData=data.frame(row.names = colnames(group_all),
                   condition=group_list1,
                   total=group_list2,
                   unique=group_list3,
                   peaks=group_list4,
                   FRiP = group_list5,
                   NSC = group_list6,
                   RSC = group_list7,
                   PBC = group_list8,
                   gender = group_list9,
                   age = group_list10,
                   PMI = group_list11)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$unique=as.numeric(colData$unique)
colData$peaks=as.numeric(colData$peaks)
colData$FRiP=as.numeric(colData$FRiP)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$PBC=as.numeric(colData$PBC)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=factor(colData$PMI)

dds=DESeqDataSetFromMatrix(countData = group_all,
                           colData = colData,
                           design = ~ condition + age + FRiP + NSC + peaks + unique)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_all[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
all_glia = group_all[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_SCZ_Control_neuron_enhancer.csv')
write.csv(new_DEG_matrix,'Differential_SCZ_Control_neuron.csv')


##########################untrt_SCZ vs. Control_1####################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/count files/")
raw_def = read.csv('untrt_SCZ_Control_1_neuron_enhancer_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
meta = factors_untrt
meta$condition_num = c(rep(1,30),rep(2,30))
group_untrt = def
miss_1 <- c()
for (i in 1:nrow(group_untrt)) {
  if(length(which(group_untrt[i,1:30]<20)) > 0.5*ncol(group_untrt[,1:30])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_untrt)) {
  if(length(which(group_untrt[i,31:60]<20)) > 0.5*ncol(group_untrt[,31:60])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_untrt = group_untrt[-miss,]
# t_group_all = t(group_untrt)
# umap_results <- umap(t_group_all,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
PCA_results = prcomp(group_untrt, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$totalreads,method = 'pearson')
cor_con_unique = cor(meta$condition_num,meta$uniquereads,method = 'pearson')
cor_con_peak = cor(meta$condition_num,meta$peaks,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_PBC = cor(meta$condition_num,meta$PBC,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$totalreads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$totalreads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$totalreads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$totalreads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$totalreads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$totalreads,method = 'pearson')

corr_uni_unique_1 = cor(PCA$PC1,meta$uniquereads,method = 'pearson')
corr_uni_unique_2 = cor(PCA$PC2,meta$uniquereads,method = 'pearson')
corr_uni_unique_3 = cor(PCA$PC3,meta$uniquereads,method = 'pearson')
corr_uni_unique_4 = cor(PCA$PC4,meta$uniquereads,method = 'pearson')
corr_uni_unique_5 = cor(PCA$PC5,meta$uniquereads,method = 'pearson')
corr_uni_unique_6 = cor(PCA$PC6,meta$uniquereads,method = 'pearson')

corr_uni_peaks_1 = cor(PCA$PC1,meta$peaks,method = 'pearson')
corr_uni_peaks_2 = cor(PCA$PC2,meta$peaks,method = 'pearson')
corr_uni_peaks_3 = cor(PCA$PC3,meta$peaks,method = 'pearson')
corr_uni_peaks_4 = cor(PCA$PC4,meta$peaks,method = 'pearson')
corr_uni_peaks_5 = cor(PCA$PC5,meta$peaks,method = 'pearson')
corr_uni_peaks_6 = cor(PCA$PC6,meta$peaks,method = 'pearson')

corr_uni_FRiP_1 = cor(PCA$PC1,meta$FRiP,method = 'pearson')
corr_uni_FRiP_2 = cor(PCA$PC2,meta$FRiP,method = 'pearson')
corr_uni_FRiP_3 = cor(PCA$PC3,meta$FRiP,method = 'pearson')
corr_uni_FRiP_4 = cor(PCA$PC4,meta$FRiP,method = 'pearson')
corr_uni_FRiP_5 = cor(PCA$PC5,meta$FRiP,method = 'pearson')
corr_uni_FRiP_6 = cor(PCA$PC6,meta$FRiP,method = 'pearson')

corr_uni_NSC_1 = cor(PCA$PC1,meta$NSC,method = 'pearson')
corr_uni_NSC_2 = cor(PCA$PC2,meta$NSC,method = 'pearson')
corr_uni_NSC_3 = cor(PCA$PC3,meta$NSC,method = 'pearson')
corr_uni_NSC_4 = cor(PCA$PC4,meta$NSC,method = 'pearson')
corr_uni_NSC_5 = cor(PCA$PC5,meta$NSC,method = 'pearson')
corr_uni_NSC_6 = cor(PCA$PC6,meta$NSC,method = 'pearson')

corr_uni_RSC_1 = cor(PCA$PC1,meta$RSC,method = 'pearson')
corr_uni_RSC_2 = cor(PCA$PC2,meta$RSC,method = 'pearson')
corr_uni_RSC_3 = cor(PCA$PC3,meta$RSC,method = 'pearson')
corr_uni_RSC_4 = cor(PCA$PC4,meta$RSC,method = 'pearson')
corr_uni_RSC_5 = cor(PCA$PC5,meta$RSC,method = 'pearson')
corr_uni_RSC_6 = cor(PCA$PC6,meta$RSC,method = 'pearson')

corr_uni_PBC_1 = cor(PCA$PC1,meta$PBC,method = 'pearson')
corr_uni_PBC_2 = cor(PCA$PC2,meta$PBC,method = 'pearson')
corr_uni_PBC_3 = cor(PCA$PC3,meta$PBC,method = 'pearson')
corr_uni_PBC_4 = cor(PCA$PC4,meta$PBC,method = 'pearson')
corr_uni_PBC_5 = cor(PCA$PC5,meta$PBC,method = 'pearson')
corr_uni_PBC_6 = cor(PCA$PC6,meta$PBC,method = 'pearson')

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

group_list1 = meta$condition
group_list2 = meta$totalreads
group_list3 = meta$uniquereads
group_list4 = meta$peaks
group_list5 = meta$FRiP
group_list6 = meta$NSC
group_list7 = meta$RSC
group_list8 = meta$PBC
group_list9 = meta$gender
group_list10 = meta$age
group_list11 = meta$PMI

colData=data.frame(row.names = colnames(group_untrt),
                   condition=group_list1,
                   total=group_list2,
                   unique=group_list3,
                   peaks=group_list4,
                   FRiP = group_list5,
                   NSC = group_list6,
                   RSC = group_list7,
                   PBC = group_list8,
                   gender = group_list9,
                   age = group_list10,
                   PMI = group_list11)

colData$condition = factor(colData$condition)
colData$total =as.numeric(colData$total)
colData$unique=as.numeric(colData$unique)
colData$peaks=as.numeric(colData$peaks)
colData$FRiP=as.numeric(colData$FRiP)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$PBC=as.numeric(colData$PBC)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=factor(colData$PMI)

dds=DESeqDataSetFromMatrix(countData = group_untrt,
                           colData = colData,
                           design = ~condition + age + FRiP + NSC + peaks + RSC + unique)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_untrt[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
all_glia = group_untrt[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_untrt_SCZ_Control_1_neuron_enhancer.csv')
write.csv(new_DEG_matrix,'Differential_untrt_SCZ_Control_1_neuron.csv')


##########################trt_SCZ vs. Control_2####################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/count files/")
raw_def = read.csv('trt_SCZ_Control_2_neuron_enhancer_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
meta = factors_trt
meta$condition_num = c(rep(1,28),rep(2,28))
group_trt = def
miss_1 <- c()
for (i in 1:nrow(group_trt)) {
  if(length(which(group_trt[i,1:28]<20)) > 0.5*ncol(group_trt[,1:28])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_trt)) {
  if(length(which(group_trt[i,29:56]<20)) > 0.5*ncol(group_trt[,29:56])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_trt = group_trt[-miss,]
# t_group_trt = t(group_trt)
# umap_results <- umap(t_group_trt,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
PCA_results = prcomp(group_trt, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$totalreads,method = 'pearson')
cor_con_unique = cor(meta$condition_num,meta$uniquereads,method = 'pearson')
cor_con_peak = cor(meta$condition_num,meta$peaks,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_PBC = cor(meta$condition_num,meta$PBC,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$totalreads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$totalreads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$totalreads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$totalreads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$totalreads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$totalreads,method = 'pearson')

corr_uni_unique_1 = cor(PCA$PC1,meta$uniquereads,method = 'pearson')
corr_uni_unique_2 = cor(PCA$PC2,meta$uniquereads,method = 'pearson')
corr_uni_unique_3 = cor(PCA$PC3,meta$uniquereads,method = 'pearson')
corr_uni_unique_4 = cor(PCA$PC4,meta$uniquereads,method = 'pearson')
corr_uni_unique_5 = cor(PCA$PC5,meta$uniquereads,method = 'pearson')
corr_uni_unique_6 = cor(PCA$PC6,meta$uniquereads,method = 'pearson')

corr_uni_peaks_1 = cor(PCA$PC1,meta$peaks,method = 'pearson')
corr_uni_peaks_2 = cor(PCA$PC2,meta$peaks,method = 'pearson')
corr_uni_peaks_3 = cor(PCA$PC3,meta$peaks,method = 'pearson')
corr_uni_peaks_4 = cor(PCA$PC4,meta$peaks,method = 'pearson')
corr_uni_peaks_5 = cor(PCA$PC5,meta$peaks,method = 'pearson')
corr_uni_peaks_6 = cor(PCA$PC6,meta$peaks,method = 'pearson')

corr_uni_FRiP_1 = cor(PCA$PC1,meta$FRiP,method = 'pearson')
corr_uni_FRiP_2 = cor(PCA$PC2,meta$FRiP,method = 'pearson')
corr_uni_FRiP_3 = cor(PCA$PC3,meta$FRiP,method = 'pearson')
corr_uni_FRiP_4 = cor(PCA$PC4,meta$FRiP,method = 'pearson')
corr_uni_FRiP_5 = cor(PCA$PC5,meta$FRiP,method = 'pearson')
corr_uni_FRiP_6 = cor(PCA$PC6,meta$FRiP,method = 'pearson')

corr_uni_NSC_1 = cor(PCA$PC1,meta$NSC,method = 'pearson')
corr_uni_NSC_2 = cor(PCA$PC2,meta$NSC,method = 'pearson')
corr_uni_NSC_3 = cor(PCA$PC3,meta$NSC,method = 'pearson')
corr_uni_NSC_4 = cor(PCA$PC4,meta$NSC,method = 'pearson')
corr_uni_NSC_5 = cor(PCA$PC5,meta$NSC,method = 'pearson')
corr_uni_NSC_6 = cor(PCA$PC6,meta$NSC,method = 'pearson')

corr_uni_RSC_1 = cor(PCA$PC1,meta$RSC,method = 'pearson')
corr_uni_RSC_2 = cor(PCA$PC2,meta$RSC,method = 'pearson')
corr_uni_RSC_3 = cor(PCA$PC3,meta$RSC,method = 'pearson')
corr_uni_RSC_4 = cor(PCA$PC4,meta$RSC,method = 'pearson')
corr_uni_RSC_5 = cor(PCA$PC5,meta$RSC,method = 'pearson')
corr_uni_RSC_6 = cor(PCA$PC6,meta$RSC,method = 'pearson')

corr_uni_PBC_1 = cor(PCA$PC1,meta$PBC,method = 'pearson')
corr_uni_PBC_2 = cor(PCA$PC2,meta$PBC,method = 'pearson')
corr_uni_PBC_3 = cor(PCA$PC3,meta$PBC,method = 'pearson')
corr_uni_PBC_4 = cor(PCA$PC4,meta$PBC,method = 'pearson')
corr_uni_PBC_5 = cor(PCA$PC5,meta$PBC,method = 'pearson')
corr_uni_PBC_6 = cor(PCA$PC6,meta$PBC,method = 'pearson')

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

group_list1 = meta$condition
group_list2 = meta$totalreads
group_list3 = meta$uniquereads
group_list4 = meta$peaks
group_list5 = meta$FRiP
group_list6 = meta$NSC
group_list7 = meta$RSC
group_list8 = meta$PBC
group_list9 = meta$gender
group_list10 = meta$age
group_list11 = meta$PMI

colData=data.frame(row.names = colnames(group_trt),
                   condition=group_list1,
                   total=group_list2,
                   unique=group_list3,
                   peaks=group_list4,
                   FRiP = group_list5,
                   NSC = group_list6,
                   RSC = group_list7,
                   PBC = group_list8,
                   gender = group_list9,
                   age = group_list10,
                   PMI = group_list11)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$unique=as.numeric(colData$unique)
colData$peaks=as.numeric(colData$peaks)
colData$FRiP=as.numeric(colData$FRiP)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$PBC=as.numeric(colData$PBC)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=factor(colData$PMI)

dds=DESeqDataSetFromMatrix(countData = group_trt,
                           colData = colData,
                           design = ~ condition + unique)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_trt[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
all_glia = group_trt[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_trt_SCZ_Control_2_neuron_enhancer.csv')
write.csv(new_DEG_matrix,'Differential_trt_SCZ_Control_2_neuron.csv')

##############################glia###########
############preparation##############
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/cofactors/")
factors = read.csv('enhancer_SCZ_Control_glia.csv')
row.names(factors) = factors[,1]
factors = factors[,-1]
factors$totalreads = round((factors$totalreads/1000000),2)
factors$uniquereads = round((factors$uniquereads/1000000),2)
factors$peaks = round((factors$peaks/1000),2)
factors$NSC = round(factors$NSC,2)
factors$RSC = round(factors$RSC,2)
factors$PBC = round(factors$PBC,2)
factors_all = factors
factors_untrt = bind_rows(factors[1:30,],factors[59:88,])
factors_trt = bind_rows(factors[31:58,],factors[89:116,])
###############################SCZ vs. Control###############################################################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/count files/")
raw_def = read.csv('SCZ_Control_glia_enhancer_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
meta = factors_all
meta$condition_num = c(rep(1,58),rep(2,58))
group_all = def
miss_1 <- c()
for (i in 1:nrow(group_all)) {
  if(length(which(group_all[i,1:58]<20)) > 0.5*ncol(group_all[,1:58])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_all)) {
  if(length(which(group_all[i,59:116]<20)) > 0.5*ncol(group_all[,59:116])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_all = group_all[-miss,]
# t_group_all = t(group_all)
# umap_results <- umap(t_group_all,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
PCA_results = prcomp(group_all, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$totalreads,method = 'pearson')
cor_con_unique = cor(meta$condition_num,meta$uniquereads,method = 'pearson')
cor_con_peak = cor(meta$condition_num,meta$peaks,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_PBC = cor(meta$condition_num,meta$PBC,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$totalreads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$totalreads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$totalreads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$totalreads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$totalreads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$totalreads,method = 'pearson')

corr_uni_unique_1 = cor(PCA$PC1,meta$uniquereads,method = 'pearson')
corr_uni_unique_2 = cor(PCA$PC2,meta$uniquereads,method = 'pearson')
corr_uni_unique_3 = cor(PCA$PC3,meta$uniquereads,method = 'pearson')
corr_uni_unique_4 = cor(PCA$PC4,meta$uniquereads,method = 'pearson')
corr_uni_unique_5 = cor(PCA$PC5,meta$uniquereads,method = 'pearson')
corr_uni_unique_6 = cor(PCA$PC6,meta$uniquereads,method = 'pearson')

corr_uni_peaks_1 = cor(PCA$PC1,meta$peaks,method = 'pearson')
corr_uni_peaks_2 = cor(PCA$PC2,meta$peaks,method = 'pearson')
corr_uni_peaks_3 = cor(PCA$PC3,meta$peaks,method = 'pearson')
corr_uni_peaks_4 = cor(PCA$PC4,meta$peaks,method = 'pearson')
corr_uni_peaks_5 = cor(PCA$PC5,meta$peaks,method = 'pearson')
corr_uni_peaks_6 = cor(PCA$PC6,meta$peaks,method = 'pearson')

corr_uni_FRiP_1 = cor(PCA$PC1,meta$FRiP,method = 'pearson')
corr_uni_FRiP_2 = cor(PCA$PC2,meta$FRiP,method = 'pearson')
corr_uni_FRiP_3 = cor(PCA$PC3,meta$FRiP,method = 'pearson')
corr_uni_FRiP_4 = cor(PCA$PC4,meta$FRiP,method = 'pearson')
corr_uni_FRiP_5 = cor(PCA$PC5,meta$FRiP,method = 'pearson')
corr_uni_FRiP_6 = cor(PCA$PC6,meta$FRiP,method = 'pearson')

corr_uni_NSC_1 = cor(PCA$PC1,meta$NSC,method = 'pearson')
corr_uni_NSC_2 = cor(PCA$PC2,meta$NSC,method = 'pearson')
corr_uni_NSC_3 = cor(PCA$PC3,meta$NSC,method = 'pearson')
corr_uni_NSC_4 = cor(PCA$PC4,meta$NSC,method = 'pearson')
corr_uni_NSC_5 = cor(PCA$PC5,meta$NSC,method = 'pearson')
corr_uni_NSC_6 = cor(PCA$PC6,meta$NSC,method = 'pearson')

corr_uni_RSC_1 = cor(PCA$PC1,meta$RSC,method = 'pearson')
corr_uni_RSC_2 = cor(PCA$PC2,meta$RSC,method = 'pearson')
corr_uni_RSC_3 = cor(PCA$PC3,meta$RSC,method = 'pearson')
corr_uni_RSC_4 = cor(PCA$PC4,meta$RSC,method = 'pearson')
corr_uni_RSC_5 = cor(PCA$PC5,meta$RSC,method = 'pearson')
corr_uni_RSC_6 = cor(PCA$PC6,meta$RSC,method = 'pearson')

corr_uni_PBC_1 = cor(PCA$PC1,meta$PBC,method = 'pearson')
corr_uni_PBC_2 = cor(PCA$PC2,meta$PBC,method = 'pearson')
corr_uni_PBC_3 = cor(PCA$PC3,meta$PBC,method = 'pearson')
corr_uni_PBC_4 = cor(PCA$PC4,meta$PBC,method = 'pearson')
corr_uni_PBC_5 = cor(PCA$PC5,meta$PBC,method = 'pearson')
corr_uni_PBC_6 = cor(PCA$PC6,meta$PBC,method = 'pearson')

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

group_list1 = meta$condition
group_list2 = meta$totalreads
group_list3 = meta$uniquereads
group_list4 = meta$peaks
group_list5 = meta$FRiP
group_list6 = meta$NSC
group_list7 = meta$RSC
group_list8 = meta$PBC
group_list9 = meta$gender
group_list10 = meta$age
group_list11 = meta$PMI

colData=data.frame(row.names = colnames(group_all),
                   condition=group_list1,
                   total=group_list2,
                   unique=group_list3,
                   peaks=group_list4,
                   FRiP = group_list5,
                   NSC = group_list6,
                   RSC = group_list7,
                   PBC = group_list8,
                   gender = group_list9,
                   age = group_list10,
                   PMI = group_list11)

colData$condition = factor(colData$condition)
colData$total =as.numeric(colData$total)
colData$unique=as.numeric(colData$unique)
colData$peaks=as.numeric(colData$peaks)
colData$FRiP=as.numeric(colData$FRiP)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$PBC=as.numeric(colData$PBC)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=factor(colData$PMI)

dds=DESeqDataSetFromMatrix(countData = group_all,
                           colData = colData,
                           design = ~ condition + NSC + RSC + unique)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_all[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
all_glia = group_all[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_SCZ_Control_glia_enhancer.csv')
write.csv(new_DEG_matrix,'Differential_SCZ_Control_glia.csv')

##########################untrt_SCZ vs. Control_1####################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/count files/")
raw_def = read.csv('untrt_SCZ_Control_1_glia_enhancer_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
meta = factors_untrt
meta$condition_num = c(rep(1,30),rep(2,30))
group_untrt = def
miss_1 <- c()
for (i in 1:nrow(group_untrt)) {
  if(length(which(group_untrt[i,1:30]<20)) > 0.5*ncol(group_untrt[,1:30])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_untrt)) {
  if(length(which(group_untrt[i,31:60]<20)) > 0.5*ncol(group_untrt[,31:60])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_untrt = group_untrt[-miss,]
# t_group_all = t(group_all)
# umap_results <- umap(t_group_all,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
PCA_results = prcomp(group_untrt, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$totalreads,method = 'pearson')
cor_con_unique = cor(meta$condition_num,meta$uniquereads,method = 'pearson')
cor_con_peak = cor(meta$condition_num,meta$peaks,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_PBC = cor(meta$condition_num,meta$PBC,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$totalreads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$totalreads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$totalreads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$totalreads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$totalreads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$totalreads,method = 'pearson')

corr_uni_unique_1 = cor(PCA$PC1,meta$uniquereads,method = 'pearson')
corr_uni_unique_2 = cor(PCA$PC2,meta$uniquereads,method = 'pearson')
corr_uni_unique_3 = cor(PCA$PC3,meta$uniquereads,method = 'pearson')
corr_uni_unique_4 = cor(PCA$PC4,meta$uniquereads,method = 'pearson')
corr_uni_unique_5 = cor(PCA$PC5,meta$uniquereads,method = 'pearson')
corr_uni_unique_6 = cor(PCA$PC6,meta$uniquereads,method = 'pearson')

corr_uni_peaks_1 = cor(PCA$PC1,meta$peaks,method = 'pearson')
corr_uni_peaks_2 = cor(PCA$PC2,meta$peaks,method = 'pearson')
corr_uni_peaks_3 = cor(PCA$PC3,meta$peaks,method = 'pearson')
corr_uni_peaks_4 = cor(PCA$PC4,meta$peaks,method = 'pearson')
corr_uni_peaks_5 = cor(PCA$PC5,meta$peaks,method = 'pearson')
corr_uni_peaks_6 = cor(PCA$PC6,meta$peaks,method = 'pearson')

corr_uni_FRiP_1 = cor(PCA$PC1,meta$FRiP,method = 'pearson')
corr_uni_FRiP_2 = cor(PCA$PC2,meta$FRiP,method = 'pearson')
corr_uni_FRiP_3 = cor(PCA$PC3,meta$FRiP,method = 'pearson')
corr_uni_FRiP_4 = cor(PCA$PC4,meta$FRiP,method = 'pearson')
corr_uni_FRiP_5 = cor(PCA$PC5,meta$FRiP,method = 'pearson')
corr_uni_FRiP_6 = cor(PCA$PC6,meta$FRiP,method = 'pearson')

corr_uni_NSC_1 = cor(PCA$PC1,meta$NSC,method = 'pearson')
corr_uni_NSC_2 = cor(PCA$PC2,meta$NSC,method = 'pearson')
corr_uni_NSC_3 = cor(PCA$PC3,meta$NSC,method = 'pearson')
corr_uni_NSC_4 = cor(PCA$PC4,meta$NSC,method = 'pearson')
corr_uni_NSC_5 = cor(PCA$PC5,meta$NSC,method = 'pearson')
corr_uni_NSC_6 = cor(PCA$PC6,meta$NSC,method = 'pearson')

corr_uni_RSC_1 = cor(PCA$PC1,meta$RSC,method = 'pearson')
corr_uni_RSC_2 = cor(PCA$PC2,meta$RSC,method = 'pearson')
corr_uni_RSC_3 = cor(PCA$PC3,meta$RSC,method = 'pearson')
corr_uni_RSC_4 = cor(PCA$PC4,meta$RSC,method = 'pearson')
corr_uni_RSC_5 = cor(PCA$PC5,meta$RSC,method = 'pearson')
corr_uni_RSC_6 = cor(PCA$PC6,meta$RSC,method = 'pearson')

corr_uni_PBC_1 = cor(PCA$PC1,meta$PBC,method = 'pearson')
corr_uni_PBC_2 = cor(PCA$PC2,meta$PBC,method = 'pearson')
corr_uni_PBC_3 = cor(PCA$PC3,meta$PBC,method = 'pearson')
corr_uni_PBC_4 = cor(PCA$PC4,meta$PBC,method = 'pearson')
corr_uni_PBC_5 = cor(PCA$PC5,meta$PBC,method = 'pearson')
corr_uni_PBC_6 = cor(PCA$PC6,meta$PBC,method = 'pearson')

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

group_list1 = meta$condition
group_list2 = meta$totalreads
group_list3 = meta$uniquereads
group_list4 = meta$peaks
group_list5 = meta$FRiP
group_list6 = meta$NSC
group_list7 = meta$RSC
group_list8 = meta$PBC
group_list9 = meta$gender
group_list10 = meta$age
group_list11 = meta$PMI

colData=data.frame(row.names = colnames(group_untrt),
                   condition=group_list1,
                   total=group_list2,
                   unique=group_list3,
                   peaks=group_list4,
                   FRiP = group_list5,
                   NSC = group_list6,
                   RSC = group_list7,
                   PBC = group_list8,
                   gender = group_list9,
                   age = group_list10,
                   PMI = group_list11)

colData$condition = factor(colData$condition)
colData$total =as.numeric(colData$total)
colData$unique=as.numeric(colData$unique)
colData$peaks=as.numeric(colData$peaks)
colData$FRiP=as.numeric(colData$FRiP)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$PBC=as.numeric(colData$PBC)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=factor(colData$PMI)

dds=DESeqDataSetFromMatrix(countData = group_untrt,
                           colData = colData,
                           design = ~ condition + age + FRiP + NSC)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_untrt[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
all_glia = group_untrt[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_untrt_SCZ_Control_1_glia_enhancer.csv')
write.csv(new_DEG_matrix,'Differential_untrt_SCZ_Control_1_glia.csv')


##########################trt_SCZ vs. Control_2####################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/count files/")
raw_def = read.csv('trt_SCZ_Control_2_glia_enhancer_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
meta = factors_trt
meta$condition_num = c(rep(1,28),rep(2,28))
group_trt = def
miss_1 <- c()
for (i in 1:nrow(group_trt)) {
  if(length(which(group_trt[i,1:28]<20)) > 0.5*ncol(group_trt[,1:28])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_trt)) {
  if(length(which(group_trt[i,29:56]<20)) > 0.5*ncol(group_trt[,29:56])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_trt = group_trt[-miss,]
# t_group_all = t(group_all)
# umap_results <- umap(t_group_all,n_neighbors=20,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
PCA_results = prcomp(group_trt, scale. = TRUE)
PCA = data.frame(PCA_results$rotation)
PCA = PCA[,1:6]

cor_con_total = cor(meta$condition_num,meta$totalreads,method = 'pearson')
cor_con_unique = cor(meta$condition_num,meta$uniquereads,method = 'pearson')
cor_con_peak = cor(meta$condition_num,meta$peaks,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_PBC = cor(meta$condition_num,meta$PBC,method = 'pearson')
cor_con_gender = cor(meta$condition_num,meta$gender,method = 'pearson')
cor_con_age = cor(meta$condition_num,meta$age,method = 'pearson')
cor_con_PMI = cor(meta$condition_num,meta$PMI,method = 'pearson')

corr_uni_total_1 = cor(PCA$PC1,meta$totalreads,method = 'pearson')
corr_uni_total_2 = cor(PCA$PC2,meta$totalreads,method = 'pearson')
corr_uni_total_3 = cor(PCA$PC3,meta$totalreads,method = 'pearson')
corr_uni_total_4 = cor(PCA$PC4,meta$totalreads,method = 'pearson')
corr_uni_total_5 = cor(PCA$PC5,meta$totalreads,method = 'pearson')
corr_uni_total_6 = cor(PCA$PC6,meta$totalreads,method = 'pearson')

corr_uni_unique_1 = cor(PCA$PC1,meta$uniquereads,method = 'pearson')
corr_uni_unique_2 = cor(PCA$PC2,meta$uniquereads,method = 'pearson')
corr_uni_unique_3 = cor(PCA$PC3,meta$uniquereads,method = 'pearson')
corr_uni_unique_4 = cor(PCA$PC4,meta$uniquereads,method = 'pearson')
corr_uni_unique_5 = cor(PCA$PC5,meta$uniquereads,method = 'pearson')
corr_uni_unique_6 = cor(PCA$PC6,meta$uniquereads,method = 'pearson')

corr_uni_peaks_1 = cor(PCA$PC1,meta$peaks,method = 'pearson')
corr_uni_peaks_2 = cor(PCA$PC2,meta$peaks,method = 'pearson')
corr_uni_peaks_3 = cor(PCA$PC3,meta$peaks,method = 'pearson')
corr_uni_peaks_4 = cor(PCA$PC4,meta$peaks,method = 'pearson')
corr_uni_peaks_5 = cor(PCA$PC5,meta$peaks,method = 'pearson')
corr_uni_peaks_6 = cor(PCA$PC6,meta$peaks,method = 'pearson')

corr_uni_FRiP_1 = cor(PCA$PC1,meta$FRiP,method = 'pearson')
corr_uni_FRiP_2 = cor(PCA$PC2,meta$FRiP,method = 'pearson')
corr_uni_FRiP_3 = cor(PCA$PC3,meta$FRiP,method = 'pearson')
corr_uni_FRiP_4 = cor(PCA$PC4,meta$FRiP,method = 'pearson')
corr_uni_FRiP_5 = cor(PCA$PC5,meta$FRiP,method = 'pearson')
corr_uni_FRiP_6 = cor(PCA$PC6,meta$FRiP,method = 'pearson')

corr_uni_NSC_1 = cor(PCA$PC1,meta$NSC,method = 'pearson')
corr_uni_NSC_2 = cor(PCA$PC2,meta$NSC,method = 'pearson')
corr_uni_NSC_3 = cor(PCA$PC3,meta$NSC,method = 'pearson')
corr_uni_NSC_4 = cor(PCA$PC4,meta$NSC,method = 'pearson')
corr_uni_NSC_5 = cor(PCA$PC5,meta$NSC,method = 'pearson')
corr_uni_NSC_6 = cor(PCA$PC6,meta$NSC,method = 'pearson')

corr_uni_RSC_1 = cor(PCA$PC1,meta$RSC,method = 'pearson')
corr_uni_RSC_2 = cor(PCA$PC2,meta$RSC,method = 'pearson')
corr_uni_RSC_3 = cor(PCA$PC3,meta$RSC,method = 'pearson')
corr_uni_RSC_4 = cor(PCA$PC4,meta$RSC,method = 'pearson')
corr_uni_RSC_5 = cor(PCA$PC5,meta$RSC,method = 'pearson')
corr_uni_RSC_6 = cor(PCA$PC6,meta$RSC,method = 'pearson')

corr_uni_PBC_1 = cor(PCA$PC1,meta$PBC,method = 'pearson')
corr_uni_PBC_2 = cor(PCA$PC2,meta$PBC,method = 'pearson')
corr_uni_PBC_3 = cor(PCA$PC3,meta$PBC,method = 'pearson')
corr_uni_PBC_4 = cor(PCA$PC4,meta$PBC,method = 'pearson')
corr_uni_PBC_5 = cor(PCA$PC5,meta$PBC,method = 'pearson')
corr_uni_PBC_6 = cor(PCA$PC6,meta$PBC,method = 'pearson')

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

group_list1 = meta$condition
group_list2 = meta$totalreads
group_list3 = meta$uniquereads
group_list4 = meta$peaks
group_list5 = meta$FRiP
group_list6 = meta$NSC
group_list7 = meta$RSC
group_list8 = meta$PBC
group_list9 = meta$gender
group_list10 = meta$age
group_list11 = meta$PMI

colData=data.frame(row.names = colnames(group_trt),
                   condition=group_list1,
                   total=group_list2,
                   unique=group_list3,
                   peaks=group_list4,
                   FRiP = group_list5,
                   NSC = group_list6,
                   RSC = group_list7,
                   PBC = group_list8,
                   gender = group_list9,
                   age = group_list10,
                   PMI = group_list11)

colData$condition = factor(colData$condition,levels = c('Control','SCZ'))
colData$total =as.numeric(colData$total)
colData$unique=as.numeric(colData$unique)
colData$peaks=as.numeric(colData$peaks)
colData$FRiP=as.numeric(colData$FRiP)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$PBC=as.numeric(colData$PBC)
colData$gender=as.numeric(colData$gender)
colData$age=as.numeric(colData$age)
colData$PMI=factor(colData$PMI)

dds=DESeqDataSetFromMatrix(countData = group_trt,
                           colData = colData,
                           design = ~ condition + age + unique)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","Control","SCZ"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.26),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_trt[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/enhancer_count/differential_peaks/")
all_glia = group_trt[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_trt_SCZ_Control_2_glia_enhancer.csv')
write.csv(new_DEG_matrix,'Differential_trt_SCZ_Control_2_glia.csv')

