rm(list = ls())
options(stringsAsFactors = F)

setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/differential analysis/ChIP-seq/promoters/diff_peaks")
def = read.csv('Differential_SCZ_Control_glia.csv')
def = def[order(def$X),]
def_data = def[,2:5]
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/differential analysis/ChIP-seq/promoters/annotated_peaks")
write.table(def_data,"Diff_SCZ_Control_glia_peaks.bed",
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "Diff_SCZ_Control_glia_peaks.bed"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")
peakAnno_df = as.data.frame(peakAnno)
data = data.frame(def[,2:5],def[,11],def[,7])
data$gene_id = 'none'
for (i in 1:nrow(data)) {
  data$gene_id[i] = peakAnno_df$SYMBOL[which(peakAnno_df$V4 == data$peak_num[i])]
}
colnames(data)[5:6] = c('adj.P','log2FC')
write.csv(data,'annotated_Diff_SCZ_Control_glia_peaks.csv',row.names = FALSE)

