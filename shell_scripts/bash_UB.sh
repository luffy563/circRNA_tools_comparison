#!/bin/bash

# UROBORUS 
path=/home/liuhongfei/tophat_fusion
for i in $(cat SRR_list.txt)
do
mkdir -p /home/liuhongfei/output/UROBORUS/$i
samtools view $path/$i/unmapped.bam > $path/$i/unmapped.sam
perl software/UROBORUS-master/bin/UROBORUS.pl -index /home/data/liuhongfei/pred_circ/bw_index/genome -gtf \
/home/data/liuhongfei/pred_circ/ref_genome/gh38.gtf -fasta /home/data/liuhongfei/pred_circ/ref_genome/genome \
$path/$i/unmapped.sam /home/liuhongfei/output/UROBORUS/$i/accepted_hits.bam
done


