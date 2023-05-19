st = ls())
options(stringsAsFactors = F)

library(fitdistrplus)
library(dplyr)
setwd('/work/cascades/bz10/VCU_human/data/HiC/new_chr/background')
chr1 = read.table('chr1.bed')
colnames(chr1)  =c('chr_loci1','start_loci1','end_loci1','chr_loci2','start_loci2','end_loci2','read','num')
head(chr1)
f=list.dirs(path = '/work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1',full.names = TRUE, recursive = FALSE)

find_promoter = function(dirpath){
  setwd(dirpath)
  Emp_P_1 = file.size('promoter_1.bed') != 0L
  Emp_P_2 = file.size('promoter_2.bed') != 0L
  if(Emp_P_1 & Emp_P_2){
    print('enter if function')
    print('both files have information')
    p_1 = read.table('promoter_1.bed')
    p_2 = read.table('promoter_2.bed')
    colnames(p_1) = c('chr_loci','start_loci','end_loci','num')
    colnames(p_2) = c('chr_loci','start_loci','end_loci','num')
    head(p_1)
    head(p_2)
    promoter_1 = chr1[which(chr1$num %in% p_1$num),]
    promoter_2 = data.frame(chr1[which(chr1$num %in% p_2$num),4:6],
                            chr1[which(chr1$num %in% p_2$num),1:3],
                            chr1[which(chr1$num %in% p_2$num),7:8])
    colnames(promoter_2)  =c('chr_loci1','start_loci1','end_loci1','chr_loci2','start_loci2','end_loci2','read','num')
    promoter = bind_rows(promoter_1,promoter_2)
    new_promoter = promoter[which(! promoter$read == 0),]
    my_data = new_promoter$read
    fit_w  <- fitdist(my_data, "weibull")
    scale = fit_w$estimate[[2]][1]
    shape = fit_w$estimate[[1]][1]
    pvalue= pweibull(my_data,shape,scale,lower.tail = F)
    FDR = p.adjust(pvalue,method = "fdr", n = length(pvalue))
    
    tmp = data.frame(my_data,pvalue,FDR)
    choose_x = tmp[which(tmp$FDR < 0.01),1]
    choose_promoter = new_promoter[which(new_promoter$read %in% choose_x),]
    choose_promoter$FDR = tmp[which(tmp$FDR < 0.01),3]
    write.table(choose_promoter,"choose_promoter.bed",row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
  } else if(Emp_P_1 & (!Emp_P_2)) {
    print('only file 1 has information')
    p_1 = read.table('promoter_1.bed')
    colnames(p_1) = c('chr_loci','start_loci','end_loci','num')
    head(p_1)
    promoter_1 = chr1[which(chr1$num %in% p_1$num),]
    new_promoter = promoter_1[which(! promoter_1$read == 0),]
    my_data = new_promoter$read
    fit_w  <- fitdist(my_data, "weibull")
    scale = fit_w$estimate[[2]][1]
    shape = fit_w$estimate[[1]][1]
    pvalue= pweibull(my_data,shape,scale,lower.tail = F)
    FDR = p.adjust(pvalue,method = "fdr", n = length(pvalue))
    tmp = data.frame(my_data,pvalue,FDR)
    choose_x = tmp[which(tmp$FDR < 0.01),1]
    choose_promoter = new_promoter[which(new_promoter$read %in% choose_x),]
    choose_promoter$FDR = tmp[which(tmp$FDR < 0.01),3]
    write.table(choose_promoter,"choose_promoter.bed",row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
  } else if(Emp_P_2 & (!Emp_P_1)){
    print('only file 2 has information')
    p_2 = read.table('promoter_2.bed')
    colnames(p_2) = c('chr_loci','start_loci','end_loci','num')
    head(p_2)
    promoter_2 = data.frame(chr1[which(chr1$num %in% p_2$num),4:6],
                            chr1[which(chr1$num %in% p_2$num),1:3],
                            chr1[which(chr1$num %in% p_2$num),7:8])
    colnames(promoter_2)  =c('chr_loci1','start_loci1','end_loci1','chr_loci2','start_loci2','end_loci2','read','num')
    new_promoter = promoter_2[which(! promoter_2$read == 0),]
    my_data = new_promoter$read
    fit_w  <- fitdist(my_data, "weibull")
    scale = fit_w$estimate[[2]][1]
    shape = fit_w$estimate[[1]][1]
    pvalue= pweibull(my_data,shape,scale,lower.tail = F)
    FDR = p.adjust(pvalue,method = "fdr", n = length(pvalue))
    tmp = data.frame(my_data,pvalue,FDR)
    choose_x = tmp[which(tmp$FDR < 0.01),1]
    choose_promoter = new_promoter[which(new_promoter$read %in% choose_x),]
    choose_promoter$FDR = tmp[which(tmp$FDR < 0.01),3]
    write.table(choose_promoter,"choose_promoter.bed",row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
  } else {
    resutls = data.frame(Emp_P_1,Emp_P_2)
    print('neither of the files has information')
    write.table(resutls,file = "no_promoter.bed",row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
    print('txt file generated')
  }
}
tmp = lapply(f, find_promoter)

