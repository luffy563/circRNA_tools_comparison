#!/bin/bash


path=software/CIRI_v2.0.6
# indexing
# bwa index -a gh38 ref_genome/genome.fa

for i in $(cat SRR_list.txt)
do
mkdir -p output/CIRI2/${i}/Alignment

bwa mem -t 20 -T 19 bwa_index/gh38.fa SRA/SRR9077/${i}.fastq.gz > output/CIRI2/${i}/Alignment/aln-se.sam
perl $path/CIRI2.pl -I output/CIRI2/${i}/Alignment/aln-se.sam -O output/CIRI2/${i}/out.ciri -F ref_genome/genome.fa -A ref_genome/gh38.gtf

done

for i in $(cat SRR_list-1.txt)
do
mkdir -p output/CIRI2/${i}/Alignment

bwa mem -t 20 -T 19 bwa_index/gh38.fa SRA/${i}/${i}_1.fastq.gz SRA/${i}/${i}_2.fastq.gz > output/CIRI2/${i}/Alignment/aln-pe.sam
perl $path/CIRI2.pl -I output/CIRI2/${i}/Alignment/aln-pe.sam -O output/CIRI2/${i}/out.ciri -F ref_genome/genome.fa -A ref_genome/gh38.gtf

done
