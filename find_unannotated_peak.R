library(ChIPseeker)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)

rm(list = ls())
options(stringsAsFactors = F)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/DifferentialGenes/untreated_or_treated/ChIP/bedfiles")
bedPeaksFile = "Differential_H3K27ac_glia_trt.bed"
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
annotated_peaks = data.frame(peakAnno_df[anno,1:3],peakAnno_df[anno,6]
                             ,peakAnno_df[anno,7],peakAnno_df[anno,17])
unanno_peaks = data.frame(peakAnno_df[!anno,1:3],peakAnno_df[!anno,6])
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/HiC/untrt_or_trt/hg38")
write.table(unanno_peaks,'unanno_peaks_H3K27ac_glia_trt.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

