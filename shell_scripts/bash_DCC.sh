#!/bin/bash


path=/home/data/liuhongfei/pred_circ
# indexing
# mkdir star_index
# STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./star_index \
# --genomeFastaFiles ref_genome/genome.fa --sjdbGTFfile ref_genome/gh38.gtf --sjdbOverhang 99
# sed -i 's/^chr//g' $path/ref_genome/ucsc_gh38_repeat.gtf
for i in $(cat SRR_list.txt)
do
mkdir -p output/DCC/${i}/Alignment
STAR --genomeDir star_index/ \
--readFilesIn $path/SRA/SRR9077/${i}.fastq.gz \
--readFilesCommand zcat \
--runThreadN 15 \
--chimSegmentMin 15 \
--chimScoreMin 15 \
--chimScoreSeparation 10 \
--chimJunctionOverhangMin 15 \
--outFilterMatchNmin 1 \
--outFilterMismatchNmax 2 \
--outReadsUnmapped Fastx \
--outSJfilterOverhangMin 15 15 15 15 \
--alignSJoverhangMin 15 \
--alignSJDBoverhangMin 15 \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 20 \
--outFilterScoreMin 1 \
--outFileNamePrefix output/DCC/${i}/Alignment/${i}
done
# -mt1 @mate1 \ # mate1 file containing the mate1 independently mapped chimeric.junction.out files
# -mt2 @mate2 \ # mate2 file containing the mate1 independently mapped chimeric.junction.out files
DCC @single_sample \ # @ is generally used to specify a file name
	  -T 20 \ 
      -D \
      -R $path/ref_genome/ucsc_gh38_repeat.gtf \ 
      -an $path/ref_genome/gh38.gtf \
      -F \
      -M \
      -Nr 1 2 \
      -fg \
      -G \
      -A $path/ref_genome/genome.fa
	  -O output/DCC

for i in $(cat SRR_list.txt)
do
mkdir -p output/DCC/${i}/Alignment
STAR --genomeDir star_index/ \
--readFilesIn $path/SRA/SRR9077/${i}.fastq.gz \
--readFilesCommand zcat \
--runThreadN 15 \
--chimSegmentMin 15 \
--chimScoreMin 15 \
--chimScoreSeparation 10 \
--chimJunctionOverhangMin 15 \
--outFilterMatchNmin 1 \
--outFilterMismatchNmax 2 \
--outReadsUnmapped Fastx \
--outSJfilterOverhangMin 15 15 15 15 \
--alignSJoverhangMin 15 \
--alignSJDBoverhangMin 15 \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 20 \
--outFilterScoreMin 1 \
--outFileNamePrefix output/DCC/${i}/Alignment/${i}
done

# -mt1 @mate1 \ # mate1 file containing the mate1 independently mapped chimeric.junction.out files
# -mt2 @mate2 \ # mate2 file containing the mate1 independently mapped chimeric.junction.out files
DCC @paired_sample \ 
	  -T 20 \ 
      -D \
      -R $path/ref_genome/ucsc_gh38_repeat.gtf \ 
      -an $path/ref_genome/gh38.gtf \
      -Pi \
      -F \
      -M \
      -Nr 1 2 \
      -fg \
      -G \
      -A $path/ref_genome/genome.fa
	  -O output/DCC
