#!/bin/bash

# find_circ
# single-end seq
for i in $(cat SRR_list.txt)
do
#mkdir output/find_circ/${i}
# Alignment
bowtie2 -p 8 --very-sensitive --score-min=C,-15,0 --mm -x bw2_index/genome -q -U SRA/SRR9077/${i}.fastq.gz | samtools view -hbuS - | samtools sort -o output/find_circ/${i}/output.bam
# Filter
samtools view -hf 4 output/find_circ/${i}/output.bam | samtools view -Sb - > output/find_circ/${i}/unmapped.bam
python software/find_circ/find_circ-1.2/unmapped2anchors.py output/find_circ/${i}/unmapped.bam | gzip > output/find_circ/${i}/anchors.fastq.gz
# Realignment
bowtie2 -p 8 --reorder --mm --score-min=C,-15,0 -q -x bw2_index/genome -U output/find_circ/${i}/anchors.fastq.gz | python software/find_circ/find_circ-1.2/find_circ.py -G ref_genome/genome.fa -p hsa_ -s output/find_circ/${i}/find_circ.sites.log > output/find_circ/${i}/find_circ.sites.bed 2> output/find_circ/${i}/find_circ.sites.reads
# Filteration
grep CIRCULAR output/find_circ/${i}/find_circ.sites.bed | grep -v chrM | awk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE > output/find_circ/${i}/circ_candidates.bed
done
# conversion
Rscript conversion.r
for i in $(cat SRR_list.txt)
do
# Annotate
CIRCexplorer2 annotate -r ref_genome/gh38_ref.txt -g ref_genome/gh38_chrID.fa -b output/find_circ/${i}/circ_candidates-1.bed -o output/find_circ/${i}/circularRNA_known.txt
done
Rscript annotation.r -s 1

# paired-end seq
for i in $(cat SRR_list-1.txt)
do
mkdir output/find_circ/${i}
# Alignment
bowtie2 -p 8 --very-sensitive --score-min=C,-15,0 --mm -x bw2_index/genome -q -1 SRA/${i}/${i}_1.fastq.gz -2 SRA/${i}/${i}_2.fastq.gz | samtools view -hbuS - | samtools sort -o output/find_circ/${i}/output.bam
# Filter
samtools view -hf 4 output/find_circ/${i}/output.bam | samtools view -Sb - > output/find_circ/${i}/unmapped.bam
python software/find_circ/find_circ-1.2/unmapped2anchors.py output/find_circ/${i}/unmapped.bam | gzip > output/find_circ/${i}/anchors.fastq.gz
# Realignment
bowtie2 -p 8 --reorder --mm --score-min=C,-15,0 -q -x bw2_index/genome -U output/find_circ/${i}/anchors.fastq.gz | python software/find_circ/find_circ-1.2/find_circ.py -G ref_genome/genome.fa -p hsa_ -s output/find_circ/${i}/find_circ.sites.log > output/find_circ/${i}/find_circ.sites.bed 2> output/find_circ/${i}/find_circ.sites.reads
# Filteration
grep CIRCULAR output/find_circ/${i}/find_circ.sites.bed | grep -v chrM | awk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE > output/find_circ/${i}/circ_candidates.bed
done

# conversion
Rscript conversion-paired.r
for i in $(cat SRR_list-1.txt)
do
# Annotate
CIRCexplorer2 annotate -r ref_genome/gh38_ref.txt -g ref_genome/gh38_chrID.fa -b output/find_circ/${i}/circ_candidates-1.bed -o output/find_circ/${i}/circularRNA_known.txt
done
Rscript annotation.r -s 0

