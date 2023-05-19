rm(list = ls())
options(stringsAsFactors = F)
library(DiffBind)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/diffbind/input_sheet")
########################################################################################################
H3K4me3 <- dba(sampleSheet = "SCZ_Control_glia_K4.csv")
H3K4me3_blacklist_remove = dba.blacklist(H3K4me3,blacklist = DBA_BLACKLIST_HG38, greylist = TRUE)
H3K4me3_consensus <- dba.peakset(H3K4me3_blacklist_remove,
                                 consensus = DBA_REPLICATE, 
                                 minOverlap = 2)
H3K4me3_consensus
H3K4me3_consensus_set <- dba(H3K4me3_consensus,
                             mask = H3K4me3_consensus$masks$Consensus, 
                             minOverlap = 1)
H3K4me3_consensus_set
consensus_peaks <- dba.peakset(H3K4me3_consensus_set, bRetrieve = TRUE)
H3K4me3_count <- dba.count(H3K4me3_blacklist_remove, 
                           summits = FALSE,
                           peaks = consensus_peaks,filter=1,
                           bScaleControl = TRUE,
                           minCount=1,
                           score=DBA_SCORE_TMM_MINUS_FULL)
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/diffbind/results")

save(H3K4me3,H3K4me3_blacklist_remove,H3K4me3_consensus,
     H3K4me3_consensus_set,consensus_peaks,H3K4me3_count,
     file = 'SCZ_Control_glia_K4.RData')

setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/diffbind/results")

H3K4me3_count_def = data.frame(H3K4me3_count[["binding"]])
H3K4me3_count[["chrmap"]]
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '1')] = 'chr1'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '5')] = 'chr10'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '6')] = 'chr11'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '7')] = 'chr12'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '8')] = 'chr13'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '9')] = 'chr14'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '14')] = 'chr15'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '15')] = 'chr16'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '16')] = 'chr17'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '20')] = 'chr18'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '21')] = 'chr19'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '22')] = 'chr2'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '24')] = 'chr20'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '25')] = 'chr21'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '26')] = 'chr22'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '32')] = 'chr3'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '33')] = 'chr4'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '35')] = 'chr5'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '37')] = 'chr6'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '38')] = 'chr7'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '39')] = 'chr8'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '40')] = 'chr9'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '72')] = 'chrX'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '73')] = 'chrY'
chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
             "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
             "chr2","chr20","chr21","chr22","chrX","chrY")
miss <- c()
for (i in 1:nrow(H3K4me3_count_def)) {
  if(!(H3K4me3_count_def$CHR[i] %in% chr_list)) miss <- append(miss,i)
}
new_H3K4me3_count_def = H3K4me3_count_def[-miss,]
new_H3K4me3_count_def[,4:ncol(new_H3K4me3_count_def)] = round(new_H3K4me3_count_def[,4:ncol(new_H3K4me3_count_def)])
new_H3K4me3_count_def = data.frame(new_H3K4me3_count_def[,1:3],
                                   rownames(new_H3K4me3_count_def),
                                   new_H3K4me3_count_def[,4:ncol(new_H3K4me3_count_def)])
colnames(new_H3K4me3_count_def)[1:4] = c('chr','start','end','peak_num')
setwd("~/Documents/Virginia Tech/Project_4_VCU_postmortem_SCZ/diffbind/count_files")
untrt_control_1 = data.frame(new_H3K4me3_count_def[,1:4],new_H3K4me3_count_def[,5:34],new_H3K4me3_count_def[,63:92])
trt_control_2 = data.frame(new_H3K4me3_count_def[,35:62],new_H3K4me3_count_def[,93:120])
write.csv(untrt_control_1,'untrt SCZ_Control 1_glia_K4.csv',row.names = FALSE)
write.csv(trt_control_2,'trt SCZ_Control 2_glia_K4.csv',row.names = FALSE)
write.csv(new_H3K4me3_count_def,'SCZ_Control_glia_K4.csv',row.names = FALSE)
write.table(new_H3K4me3_count_def[,1:4],'SCZ_Control_glia_K4.bed',sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)




