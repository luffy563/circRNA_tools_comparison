#!/bin/bash

# segemehl
index_path=/home/data1/liuhongfei/pred_circ/
mkdir -p se_index/
# indexing
# bash $path/segemehl.x -x se_index/gh38.idx -d ref_genome/genome.fa
# single-end seq
for i in $(cat SRR_list.txt)
do
mkdir -p output/segemehl/${i}
# Alignment
segemehl.x -t 20 -i se_index/gh38.idx -d ref_genome/genome.fa -q SRA/SRR9077/${i}/${i}.fastq -S -o output/segemehl/${i}
done
# conversion
# Rscript conversion.r
# for i in $(cat SRR_list.txt)
# do
# # Annotate
# CIRCexplorer2 annotate -r ref_genome/gh38_ref.txt -g ref_genome/gh38_chrID.fa -b output/find_circ/${i}/circ_candidates-1.bed -o output/find_circ/${i}/circularRNA_known.txt
# done
# Rscript annotation.r -s 1

# paired-end seq
for i in $(cat SRR_list-1.txt)
do
mkdir -p output/segemehl/${i}
# Alignment
segemehl.x -t 20 -i se_index/gh38.idx -d ref_genome/genome.fa -q SRA/${i}/${i}_1.fastq -p SRA/${i}/${i}_2.fastq -S -o output/segemehl/${i}
done

# conversion
# Rscript conversion-paired.r
# for i in $(cat SRR_list-1.txt)
# do
# # Annotate
# CIRCexplorer2 annotate -r ref_genome/gh38_ref.txt -g ref_genome/gh38_chrID.fa -b output/find_circ/${i}/circ_candidates-1.bed -o output/find_circ/${i}/circularRNA_known.txt
# done
# Rscript annotation.r -s 0

