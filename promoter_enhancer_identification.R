################################pre##########################
library(ChIPseeker)
#BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
rm(list = ls())
options(stringsAsFactors = F)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
##############H3K4me3###############
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/diffbind/count_files")
bedPeaksFile = "SCZ_Control_neuron_K4.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[anno,]

H3K4me3_peak = read.table('untrt_trt_SCZ_Control_glia_K4.bed')
H3K4me3_old = read.csv('trt SCZ_Control 2_glia_K4.csv')
H3K4me3 = data.frame(H3K4me3_peak,H3K4me3_old)
colnames(H3K4me3_peak) = c('chr','start','end','peak_num')
colnames(H3K4me3)[1:4] = c('chr','start','end','peak_num')
promoter_peak = H3K4me3_peak[which(H3K4me3_peak$peak_num %in% new_peakAnno_df$V4),]
promoter_count = H3K4me3[which(H3K4me3$peak_num %in% new_peakAnno_df$V4),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/differential analysis/ChIP-seq/promoters/count files")
write.csv(promoter_count,'promoter_trt SCZ_Control 2_glia.csv',row.names = FALSE)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/differential analysis/ChIP-seq/promoters/peak_files")
write.table(promoter_peak,'promoter_untrt_trt_SCZ_Control_glia.bed',
            row.names = FALSE, sep="\t", quote = FALSE)
######################H3K27ac####################
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/diffbind/count_files")
bedPeaksFile = "untrt_trt_SCZ_Control_glia_K27.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]

H3K4me3_peak = read.table('untrt_trt_SCZ_Control_glia_K27.bed')
H3K4me3_old = read.csv('trt SCZ_Control 2_glia_K27.csv')
H3K4me3 = data.frame(H3K4me3_peak,H3K4me3_old)
colnames(H3K4me3_peak) = c('chr','start','end','peak_num')
colnames(H3K4me3)[1:4] = c('chr','start','end','peak_num')
enhancer_peak = H3K4me3_peak[which(H3K4me3_peak$peak_num %in% new_peakAnno_df$V4),]
enhancer_count = H3K4me3[which(H3K4me3$peak_num %in% new_peakAnno_df$V4),]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/differential analysis/ChIP-seq/enhancers/count files")
write.csv(enhancer_count,'enhancer_trt SCZ_Control 2_glia.csv',row.names = FALSE)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/differential analysis/ChIP-seq/enhancers/peak_files")
write.table(enhancer_peak,'enhancer_untrt_trt_SCZ_Control_glia.bed',
            row.names = FALSE, sep="\t", quote = FALSE)

