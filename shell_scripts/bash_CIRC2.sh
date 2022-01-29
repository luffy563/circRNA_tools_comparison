#!/bin/bash

# CIRCexplorer2
path=/home/data/liuhongfei/pred_circ
# Alignment (single-end seq)
for i in $(cat SRR_list.txt)
do

mkdir -p output/tophat/${i}
mkdir -p output/tophat_fusion/${i}
tophat2 -a 6 --microexon-search -m 2 -p 20 -G $path/ref_genome/gh38.gtf -o output/tophat/${i} $path/bw2_index/genome $path/SRA/SRR9077/${i}.fastq.gz
bamToFastq -i output/tophat/${i}/unmapped.bam -fq output/tophat/${i}/unmapped.fastq
tophat2 -o output/tophat_fusion/${i} -p 20 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $path/bw_index/genome output/tophat/${i}/unmapped.fastq
# nohup CIRCexplorer2 align -G ref_genome/gh38_ref.gtf -i bw_index/genome -j bw2_index/genome -f SRA/SRR9077/${i}.fastq.gz -o output/tophat/${i}

# Parse
CIRCexplorer2 parse -t TopHat-Fusion output/tophat_fusion/${i}/accepted_hits.bam -b output/tophat_fusion/${i}/circ_candidates.bed > output/tophat_fusion/${i}/tophat_fusion.log
done
# Conversion to chrID
Rscript conversion.r -a CIRC2 -s 1
for i in $(cat SRR_list.txt)
do
# Annotate
CIRCexplorer2 annotate -r $path/ref_genome/gh38_ref.txt -g $path/ref_genome/gh38_chrID.fa -b output/tophat_fusion/${i}/circ_candidates_convert.bed -o output/tophat_fusion/${i}/circularRNA_known.txt
done
Rscript annotation.r -a CIRC2 -s 1


# CIRCexplorer2
# Alignment (paired-end seq)
for i in $(cat SRR_list-1.txt)
do
mkdir output/tophat/${i}
mkdir output/tophat_fusion/${i}
tophat2 -o tophat_fusion -p 20 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $path/bw_index/genome $path/SRA/${i}/${i}_1.fastq.gz $path/SRA/${i}/${i}_2.fastq.gz


# Parse
CIRCexplorer2 parse --pe -t TopHat-Fusion output/tophat_fusion/${i}/accepted_hits.bam > output/tophat_fusion/${i}/tophat_fusion.log
done
# Conversion to chrID
Rscript conversion.r -a CIRC2 -s 0
for i in $(cat SRR_list-1.txt)
do
# Annotate
CIRCexplorer2 annotate -r $path/ref_genome/gh38_ref.txt -g $path/ref_genome/gh38_chrID.fa -b output/tophat_fusion/${i}/back_spliced_junction.bed -o output/tophat_fusion/${i}/circularRNA_known.txt
done

Rscript annotation.r -a CIRC2 -s 0
