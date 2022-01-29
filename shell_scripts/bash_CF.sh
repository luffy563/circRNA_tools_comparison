#!/bin/bash


path=software/circRNA_finder-master
# indexing
mkdir star_index
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./star_index \
--genomeFastaFiles ref_genome/genome.fa --sjdbGTFfile ref_genome/gh38.gtf --sjdbOverhang 99

for i in $(cat SRR_list.txt)
do
mkdir -p output/CircRNAfinder/${i}/Alignment
STAR --genomeDir star_index/ \
--readFilesCommand gunzip -c \
--readFilesIn SRA/SRR9077/${i}.fastq.gz \
--runThreadN 20 \
--chimSegmentMin 20 \
--chimScoreMin 1 \
--alignIntronMax 500000 \
--outFilterMismatchNmax 4 \
--alignTranscriptsPerReadNmax 100000 \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--chimOutType SeparateSAMold \
--outFilterMultimapNmax 2 \
--outFileNamePrefix output/CircRNAfinder/${i}/Alignment/${i}

perl $path/postProcessStarAlignment.pl --starDir output/CircRNAfinder/${i} --minLen 100 --outDir output/CircRNAfinder/${i}
done