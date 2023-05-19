rm(list = ls())
options(stringsAsFactors = F)
candidate_genes = c('AKT1','APOE','BDNF','CHRNA7','COMT','DAO','DAOA','DISC1','DRD2','DRD3',
                    'DRD4','DTNBP1','GRM3','HTR2A','KCNN3','MTHFR','NOTCH4','NRG1','PPP3CC',
                    'PRODH','RGS4','SLC6A3','SLC6A4','TNF','ZDHHC8','CTNND2','NRXN3','PTPRR',
                    'SLIT2','ANK3','APP','PTK2','FYN','ROBO1','PLPPR4','GRIK2','CDK5R1','GRIA2',
                    'GRID2','GRIN3A','GRIN2B','GRIN2A','GRIA3','GRIA4','GRIK2','TIAM1','GRIA1',
                    'GRIN1','CLN3','GRID1','CPEB4','PTK2B','AKT1','APOE','BDNF','CHRNA7','COMT',
                    'DAO','DAOA','DISC1','DRD2','DRD3','DRD4','DTNBP1','GRM3','HTR2A','KCNN3',
                    'MTHFR','NOTCH4','NRG1','PPP3CC','PRODH','RGS4','SLC6A3','SLC6A4','TNF',
                    'ZDHHC8','MEF2C','MUSK','PPP1R9A','NRXN1','TPBG','LRRTM3','ROBO2','DAG1',
                    'FA2H','DICER1','NTRK2','SH3TC2','NDRG1','ARHGEF10','POU3F2','SOD1','ILK',
                    'CDK5','AKT1','APOE','BDNF','CHRNA7','COMT','DAO','DAOA','DISC1','DRD2',
                    'DRD3','DRD4','DTNBP1','GRM3','HTR2A','KCNN3','MTHFR','NOTCH4','NRG1',
                    'PPP3CC','PRODH','RGS4','SLC6A3','SLC6A4','TNF','ZDHHC8','CPLX2','SNCA',
                    'SLC30A1','NLGN1','P2RY1','SYT1','CNR1','CACNB2','PRKCB','DTNBP1','AKT1',
                    'APOE','BDNF','CHRNA7','COMT','DAO','DAOA','DISC1','DRD2','DRD3','DRD4',
                    'DTNBP1','GRM3','HTR2A','KCNN3','MTHFR','NOTCH4','NRG1','PPP3CC','PRODH',
                    'RGS4','SLC6A3','SLC6A4','TNF','ZDHHC8','EBF1','CHD9','PCK1','AKT1','APOE',
                    'BDNF','CHRNA7','COMT','DAO','DAOA','DISC1','DRD2','DRD3','DRD4','DTNBP1',
                    'GRM3','HTR2A','KCNN3','MTHFR','NOTCH4','NRG1','PPP3CC','PRODH'
                    ,'RGS4','SLC6A3','SLC6A4','TNF','ZDHHC8','EBF1','CHD9','PCK1','DRD5','TCEB2'
                    ,'CHRNA2','CHRNA6','FGF17','SOX9','GSN','PLP1','CD9','BOK','ASPA','SOX8'
                    ,'OLIG2','OLIG1','KCNJ10','CNP','SOX13','MYRF','GRM5','GRM7','GRM8','AKT1','APOE','BDNF','CHRNA7','COMT','DAO','DAOA','DISC1','DRD2','DRD3'
                    ,'DRD4','DTNBP1','GRM3','HTR2A','KCNN3','MTHFR','NOTCH4','NRG1','PPP3CC','PRODH'
                    ,'RGS4','SLC6A3','SLC6A4','TNF','ZDHHC8','EGR3','FZD3','GRIN2A','NTNG2','DVL2'
                    ,'DLG4', 'MACF1','FZD4')
tmp = unique(candidate_genes)

#SZTR_scRNA = read.csv('SZTR_excitatory_neuron_DEGs.csv')
#SZTR_scRNA_DEGs = unique(SZTR_scRNA$Genes)
#neuron_scRNA = read.csv('neuron_DEGs.csv')
#neuron_scRNA_DEGs = unique(neuron_scRNA$Gene)
glia_scRNA =read.csv('glia_DEGs.csv')
glia_scRNA_DEGs = unique(glia_scRNA$Genes)
#DEGs_mix = read.csv('differential_expressed_genes.csv')
#DEGs_mix = DEGs_mix[which(DEGs_mix$SCZ.fdr < 0.05),]

bulk_neuron = read.csv('DEG_glia_0.05.csv')

overlap = intersect(glia_scRNA_DEGs,bulk_neuron$X)

important_genes = intersect(overlap,tmp)

new_important_genes = intersect(tmp,bulk_neuron$X)
