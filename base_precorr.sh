#!/bin/bash

Group=GROUPGOESHERE
subgroup=SUBgroupGOESHERE
Histone=HISTONEGOESHERE

cd /work/cascades/bz10/VCU_mouse/data/$Group/$Histone/$subgroup

cd Raw_Data/
gunzip *.gz
trim_galore *.fastq

cd ..

FILES=$PWD/Raw_Data/*.fq
FQ=.fq
SAM=.sam
LOG=.log
BAM=.bam
SORT=_sort.bam
UNI=_unique.bam
PREBED=_pre.bed
BED=.bed
EXT=_extend.bed
EXPWIN=_extend_promotor_win.bed
EXG4000=_extend_geno_4000.bed
EXPCOR=_extend_promotor_nor.bed
EXGCOR=_extend_geno_4000_nor.bed
G100=_geno_win_100.bed
SORTG100=_geno_win_100_sort.bed
BedG=.bedGraph
NEBw=.bw

mkdir Aligned_BAM
mkdir Aligned_SAM
mkdir SORT_BAM
mkdir UNIQUE_BAM
mkdir BED
mkdir EXTEND_BED
mkdir EXT_PROMOTOR_WIN
mkdir EXT_GENO_WIN
mkdir Correlation
mkdir MACS
mkdir GenoWin100Sort
mkdir BedGraph
mkdir Nor_Ext_BW

input_length=$(wc -l < /work/cascades/bz10/input/VCU_mouse_brain_neuron/$Group/input.bed )

for fn in $FILES
do
echo `basename "$fn"`
f=`basename "${fn%.*}"`
echo $f

bowtie2 -p 16 -x /home/bz10/Data/index_bw2/mm10 -U $PWD/Raw_Data/$f$FQ -S $PWD/Aligned_SAM/$f$SAM 2>$PWD/Aligned_SAM/$f$LOG

samtools view -bq 10 $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$BAM

samtools sort $PWD/Aligned_BAM/$f$BAM -o $PWD/SORT_BAM/$f$SORT

samtools index $PWD/SORT_BAM/$f$SORT

samtools rmdup -s $PWD/SORT_BAM/$f$SORT $PWD/UNIQUE_BAM/$f$UNI

bedtools bamtobed -i $PWD/UNIQUE_BAM/$f$UNI > $PWD/BED/$f$PREBED

bedtools subtract -a $PWD/BED/$f$PREBED -b /home/bz10/blacklist/mm10.blacklist.bed > $PWD/BED/$f$BED

bedtools slop -i $PWD/BED/$f$BED -g /home/bz10/Data/Promotor/mm10/mm10genome.bed -b 100 > $PWD/EXTEND_BED/$f$EXT

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/bz10/Data/Promotor/mm10/mm10Promotor.bed > $PWD/EXT_PROMOTOR_WIN/$f$EXPWIN

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/bz10/Data/Promotor/mm10/mm10genome_win_4000.bed > $PWD/EXT_GENO_WIN/$f$EXG4000

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/bz10/Data/Promotor/mm10/mm10genome_win_100.bed > $PWD/EXT_GENO_WIN/$f$G100

sort -k1,1 -k2,2g -u -o $PWD/GenoWin100Sort/$f$SORTG100 $PWD/EXT_GENO_WIN/$f$G100

ChIP_length=$(wc -l < $PWD/BED/$f$BED)

paste $PWD/EXT_PROMOTOR_WIN/$f$EXPWIN /work/cascades/bz10/input/VCU_mouse_brain_neuron/$Group/input_extend_promotor_win.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXPCOR

paste $PWD/EXT_GENO_WIN/$f$EXG4000 /work/cascades/bz10/input/VCU_mouse_brain_neuron/$Group/input_extend_genome_4000.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXGCOR

macs2 callpeak -t $PWD/BED/$f$BED -c /work/cascades/bz10/input/VCU_mouse_brain_neuron/$Group/input.bed -f BED -g mm -n $f -q 0.05 --outdir $PWD/MACS

paste $PWD/EXT_GENO_WIN/$f$G100 /work/cascades/bz10/input/VCU_mouse_brain_neuron/$Group/input_extend_genome_100.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/BedGraph/$f$BedG

bedSort $PWD/BedGraph/$f$BedG $PWD/BedGraph/$f$BedG

bedGraphToBigWig $PWD/BedGraph/$f$BedG /home/bz10/Data/Promotor/mm10/mm10genome.bed $PWD/Nor_Ext_BW/$f$NEBw

done

/work/cascades/bz10/VCU_mouse/script/Summary.sh .
