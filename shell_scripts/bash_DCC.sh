#!/bin/bash

# DCC
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts
config_dir=$working_dir/config

## Indexing
# mkdir -p $working_dir/star_index
# STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $working_dir/star_index \
# --genomeFastaFiles $working_dir/ref_genome/genome.fa --sjdbGTFfile $working_dir/ref_genome/gh38.gtf --sjdbOverhang 99

## Convert chr* to *
# sed -i 's/^chr//g' $working_dir/ref_genome/ucsc_gh38_repeat.gtf

# Alignment (paired end)
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p $working_dir/output/DCC/${i}/Alignment

STAR --genomeDir $working_dir/star_index/ \
--readFilesIn $working_dir/SRA/${i}/${i}_1.fastq $working_dir/SRA/${i}/${i}_2.fastq \
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
--outFileNamePrefix $working_dir/output/DCC/${i}/Alignment/${i}
echo $working_dir/output/DCC/${i}/Alignment/${i}Chimeric.out.junction >>$config_dir/paired_sample
done

DCC @$config_dir/paired_sample \ 
	  -T 20 \ 
      -D \
      -R $working_dir/ref_genome/ucsc_gh38_repeat.gtf \ 
      -an $working_dir/ref_genome/gh38.gtf \
      -Pi \
      -F \
      -M \
      -Nr 1 2 \
      -fg \
      -G \
      -A $working_dir/ref_genome/genome.fa
	  -O $working_dir/output/DCC

# Conversion
Rscript $working_dir/conversion.r -w $working_dir -a DCC -i $working_dir/SRR_list.txt

# Annotate
for i in $(cat $working_dir/SRR_list.txt)
do

CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b output/DCC/${i}/circ_candidates_convert.bed -o $working_dir/output/DCC/${i}/circularRNA_known.txt
done

# Merged
Rscript $working_dir/annotation.r -w $working_dir -a DCC -i $working_dir/SRR_list.txt