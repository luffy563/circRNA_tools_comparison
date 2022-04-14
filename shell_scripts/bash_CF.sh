#!/bin/bash

# CircRNAfinder
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
path=$working_dir/software/circRNA_finder-master
script_dir=$working_dir/scripts

# Indexing
mkdir -p $working_dir/star_index
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $working_dir/star_index \
--genomeFastaFiles $working_dir/ref_genome/genome.fa --sjdbGTFfile $working_dir/ref_genome/gh38.gtf --sjdbOverhang 99

cat $working_dir/readsLen.txt | while read i len
do
mkdir -p $working_dir/output/CircRNAfinder/${i}/Alignment
STAR --genomeDir $working_dir/star_index/ \
# --readFilesCommand gunzip -c \
--readFilesIn $working_dir/SRA/${i}/${i}_1.fastq $working_dir/SRA/${i}/${i}_2.fastq \
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

perl $path/postProcessStarAlignment.pl --starDir $working_dir/output/CircRNAfinder/${i} --minLen $len --outDir $working_dir/output/CircRNAfinder/${i}
done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a CF -i $working_dir/SRR_list.txt

# Annotate
for i in $(cat $working_dir/SRR_list.txt)
do
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/CircRNAfinder/${i}/circ_candidates_convert.bed -o $working_dir/output/CircRNAfinder/${i}/circularRNA_known.txt
done

# Merged
Rscript $script_dir/annotation.r -w $working_dir -a CF -i $working_dir/SRR_list.txt