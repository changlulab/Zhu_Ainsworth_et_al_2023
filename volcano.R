####################prep###############
rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
library(ggrepel)
library(dplyr)
candidate_genes = c('AKT1','APOE',"BDNF","CHRNA7","COMT","DAO","DAOA","DISC1","DRD2","DRD3",
                    "DRD4","DTNBP1","GRM3","HTR2A","KCNN3","MTHFR","NOTCH4","NRG1","PPP3CC,PRODH",
                    "RGS4","SLC6A3","SLC6A4","TNF","ZDHHC8","CTNND2","NRXN3","PTPRR","SLIT2","ANK3",
                    "APP","PTK2","FYN","ROBO1","PLPPR4","GRIK2","CDK5R1","GRIA2","GRID2","GRIN3A",
                    "GRIN2B","GRIN2A","GRIA3","GRIA4","TIAM1","GRIA1","GRIN1","CLN3","GRID1","CPEB4",
                    "PTK2B","MEF2C","MUSK","PPP1R9A","NRXN1","TPBG","LRRTM3","ROBO2","DAG1","FA2H",
                    "DICER1","NTRK2","SH3TC2","NDRG1","ARHGEF10","POU3F2","SOD1","ILK","CDK5","CPLX2",
                    "SNCA","SLC30A1","NLGN1","P2RY1","SYT1","CNR1","CACNB2","PRKCB","EBF1","CHD9",
                    "PCK1","DRD5","TCEB2","CHRNA2","CHRNA6","FGF17","SOX9","GSN","PLP1","CD9",
                    "BOK","ASPA","SOX8","OLIG2","OLIG1","KCNJ10","CNP","SOX13","MYRF","GRM5",
                    "GRM7","GRM8","EGR3","NTNG2","DVL2","DLG4","MACF1","FZD4",'SHANK3','VIPR2',
                    'SLC1A1','RBM12',"APOL2","SCZD12","CHI3L1","DISC2","SYN2","SCZD1","SCZD3","SCZD5",
                    "SCZD6","SCZD11","SCZD2","SCZD7","SCZD8","RTN4R","APOL4",'G72','PRODH2','ZDHHC8',
                    'TAAR6','NOS1AP','TPH1','AHI1','RELN','TRKA','ZNF804','TCF','NRGRN','BCL9','NRGN',
                    'TCF4','ARVCF','CNNM2','CACNA1C')
tmp = unique(candidate_genes)
##############################################H3K27ac_neuron##########################################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/promoter_count/differential_peaks/")
peakAnno_df = read.csv('All_SCZ_Control_glia_promoter.csv')
H3K4me3_neuron = data.frame(peakAnno_df[,1:7])
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/ChIP-seq_differential/promoter_count/annotated_peaks")
annotated_peaks = read.csv('anno_SCZ_Control_glia_promoter.csv')
H3K4me3_neuron$SYMBOL = 'name'
for (i in 1:nrow(H3K4me3_neuron)) {
  if (H3K4me3_neuron$X[i] %in% annotated_peaks$peak_num) {
    H3K4me3_neuron$SYMBOL[i] = annotated_peaks$gene_id[which(annotated_peaks$peak_num == H3K4me3_neuron$X[i])]
  }
}
H3K4me3_neuron$diffexper <- "NO"
H3K4me3_neuron$diffexper[H3K4me3_neuron$log2FoldChange> 0.26 & H3K4me3_neuron$padj<0.05] <- "up"
H3K4me3_neuron$diffexper[H3K4me3_neuron$log2FoldChange< -0.26 & H3K4me3_neuron$padj<0.05] <- "down"
H3K4me3_neuron$genelable = 'no'

b= H3K4me3_neuron[H3K4me3_neuron$diffexper == 'up' | H3K4me3_neuron$diffexper == 'down',]
b_candidate = b[b$SYMBOL %in% tmp,]
b_candidate$row = row.names(b_candidate)
b_candidate = b_candidate %>%
  group_by(SYMBOL) %>%
  filter(padj == min(padj))
rownames(b_candidate) = b_candidate$row
a_candidate = row.names(b_candidate)
for (i in a_candidate) {
  n=as.numeric(i)
  H3K4me3_neuron$genelable[n] <- H3K4me3_neuron$SYMBOL[n]
}

point_color = c("#CCCCCC","#FF99CC","#FF99CC")
names(point_color) = c("NO","up","down")
colnames(H3K4me3_neuron)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/Cell reports_review/figures for revision/volcano_plot")
write.csv(H3K4me3_neuron,'SCZ_Control_glia_promoter.csv')

pdf('volcano_glia_promoter.pdf', width = 10, height = 8)
ggplot(data = H3K4me3_neuron,aes(x=log2FoldChange, y= -log10(padj), color=diffexper)) +
  geom_point() +
  scale_color_manual(values = point_color) +
  geom_hline(yintercept = -log10(0.05),col="black",linetype = "dashed") +
  geom_point(data = subset(H3K4me3_neuron,genelable != 'no'), color='black',size=2) +
  geom_text_repel(data = H3K4me3_neuron[a_candidate,],
                  aes(x=log2FoldChange, y= -log10(padj),label=genelable),
                  color='black',size=6,min.segment.length = 0,point.padding = 0,box.padding = 0.4) +
  xlim(-2,2) + ylim(0,10)
  
dev.off()
#H3K27ac: #66CC00
#H3K4me3: #FF99CC
#DEGs: #FFCC66
#########################################DEG_glia###################################################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/RNA_DEseq/UCSD_gtf/neuron/")
DEG_glia = read.csv('All_SCZ_Control_neuron.csv')
rownames(DEG_glia) = DEG_glia[,1]
DEG_glia = DEG_glia[,-1]
DEG_glia = DEG_glia[,1:6]
DEG_glia$gene = row.names(DEG_glia)
row.names(DEG_glia) = as.numeric(c(1:nrow(DEG_glia)))
DEG_glia$diffexper <- "NO"
DEG_glia$diffexper[DEG_glia$log2FoldChange > 0.26 & DEG_glia$padj<0.05] <- "up"
DEG_glia$diffexper[DEG_glia$log2FoldChange < -0.26 & DEG_glia$padj<0.05] <- "down"
DEG_glia$genelable = 'gene'

b= DEG_glia[DEG_glia$diffexper == 'up' | DEG_glia$diffexper == 'down',]
#b_candidate = b[b$gene %in% candidate_genes,]
b_candidate = b[b$gene %in% tmp,]
b_candidate$row = row.names(b_candidate)
b_candidate = b_candidate %>%
  group_by(gene) %>%
  filter(padj == min(padj))
rownames(b_candidate) = b_candidate$row
a_candidate = as.numeric(row.names(b_candidate))

DEG_glia$genelable[a_candidate] <- DEG_glia$gene[a_candidate]
point_color = c("#CCCCCC","#FFCC66","#FFCC66")
names(point_color) = c("NO","up","down")
#colnames(DEG_glia)
#class(DEG_glia$Conc_untreated.SCZ)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/Cell reports_review/figures for revision/volcano_plot")
write.csv(DEG_glia,'SCZ_Control_neuron_DEG.csv')

pdf('volcano_neuron_DEG.pdf', width = 10, height = 8)
ggplot(data = DEG_glia,aes(x=log2FoldChange, y= -log10(padj), color=diffexper)) +
  geom_point() +
  scale_color_manual(values = point_color) +
  geom_hline(yintercept = -log10(0.05),col="black",linetype = "dashed") +
  geom_point(data = DEG_glia[a_candidate,], color='black',size=2) +
  geom_text_repel(data = DEG_glia[a_candidate,],
                  aes(x=log2FoldChange, y= -log10(padj),label=genelable),
                  color='black',size=6,min.segment.length = 0,point.padding = 0,box.padding = 0.4) +
  ylim(0,10) + xlim(-2.5,2.5)
dev.off()





