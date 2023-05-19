 #!/bin/bash

cd /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1

mkdir bedfiles_1
mkdir bedfiles_2

cd chr1_group1
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group2
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group3
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group4
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group5
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group6
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group7
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group8
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ../chr1_group9
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/choose_promoter.bed ../bedfiles_1/$id\_promoter.bed);done
ls -d chr1_*/ | cut -d '/' -f 1| while read id;do(cp ./$id/new_choose_promoter.bed ../bedfiles_2/$id\_promoter.bed);done

cd ..
mkdir raw_loops
mkdir raw_loops/method_1
mkdir raw_loops/method_2

cd bedfiles_1
cat *.bed | grep -v 'chr_loci1'  > ../raw_loops/method_1/Loops.bed

cd ../bedfiles_2
cat *.bed | grep -v 'chr_loci1'  > ../raw_loops/method_2/Loops.bed

cd ../raw_loops

cp /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/promoter/chr1.bed /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/raw_loops/promoter.bed

Rscript /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/refine_loops.r /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/raw_loops/method_1

cd method_1
bedtools intersect -wa -wb -f 0.50 -a ../promoter.bed -b loops_part1.bed  > overlap_loops_part1.bed
bedtools intersect -wa -wb -f 0.50 -a ../promoter.bed -b loops_part2.bed  > overlap_loops_part2.bed

cd ../../ && mkdir final_loops

mkdir final_loops/method_1
mkdir final_loops/method_2

Rscript /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/final_loops.r /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/raw_loops/method_1 /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/final_loops/method_1

cd raw_loops

Rscript /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/refine_loops.r /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/raw_loops/method_2

cd method_2
bedtools intersect -wa -wb -f 0.50 -a ../promoter.bed -b loops_part1.bed  > overlap_loops_part1.bed
bedtools intersect -wa -wb -f 0.50 -a ../promoter.bed -b loops_part2.bed  > overlap_loops_part2.bed

Rscript /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/final_loops.r /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/raw_loops/method_2 /work/cascades/bz10/VCU_human/data/HiC/promoter_loops/chr1/final_loops/method_2
