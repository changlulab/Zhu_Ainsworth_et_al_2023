#!/bin/bash

echo "choose the group you want to process and press [ENTER] :"
ls /work/cascades/bz10/VCU_mouse/data
read GROUP
cd /work/cascades/bz10/VCU_mouse/data/$GROUP 2> /dev/null

##############################################################################
until [ "$PWD" == "/work/cascades/bz10/VCU_mouse/data/$GROUP" ]; do
#the space in the bracket (at the begining and at the end) is very important##

echo "Please choose a real group"
ls /work/cascades/bz10/VCU_mouse/data
read GROUP
cd /work/cascades/bz10/VCU_mouse/data/$GROUP 2> /dev/null;
done

echo "choose the hisotne you want to process and press [ENTER] :"
ls /work/cascades/bz10/VCU_mouse/data/$GROUP
read HISTONE
cd /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE 2> /dev/null
until [ "$PWD" == "/work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE" ]; do
echo "Please choose a real histone mark"
ls /work/cascades/bz10/VCU_mouse/data/$GROUP
read HISTONE
cd /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE 2> /dev/null;
done

echo "choose the subgroup you want to process and press [ENTER] :"
ls /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE
read SUBGROUP
cd /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP 2> /dev/null
until [ "$PWD" == "/work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP" ]; do
echo "Please choose a real subgroup"
ls /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE
read SUBGROUP
cd /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP 2> /dev/null;
done

cd /work/cascades/bz10/Raw_Data/MIA/ChIP_seq
sed "s/GROUPGOESHERE/"$GROUP"/" ./base_precorr.sh > /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP/precorr_1.sh

sed "s/SUBgroupGOESHERE/"$SUBGROUP"/" /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP/precorr_1.sh > /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP/precorr_2.sh

sed "s/HISTONEGOESHERE/"$HISTONE"/" /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP/precorr_2.sh > /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP/precorr.sh

cd /work/cascades/bz10/VCU_mouse/data/$GROUP/$HISTONE/$SUBGROUP

rm ./precorr_1.sh ./precorr_2.sh
chmod u+x ./precorr.sh
sbatch -A chipseq --mem-per-cpu=20G -t 3-00:00:00 -p normal_q ./precorr.sh
