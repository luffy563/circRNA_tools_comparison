#!/bin/bash

# CIRI2
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts
path=$working_dir/software/CIRI_v2.0.6

# Indexing
bwa index -a gh38 $working_dir/ref_genome/genome.fa

# Alignment (paired-end seq)
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p $working_dir/output/CIRI2/${i}/Alignment

bwa mem -t 20 -T 19 $working_dir/bwa_index/gh38.fa $working_dir/SRA/${i}/${i}_1.fastq.gz $working_dir/SRA/${i}/${i}_2.fastq.gz > $working_dir/output/CIRI2/${i}/Alignment/aln-se.sam
perl $path/CIRI2.pl -I $working_dir/output/CIRI2/${i}/Alignment/aln-se.sam -O $working_dir/output/CIRI2/${i}/out.ciri -F $working_dir/ref_genome/genome.fa -A $working_dir/ref_genome/gh38.gtf

done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a CIRI2 -i $working_dir/SRR_list.txt

# Annotate
for i in $(cat $working_dir/SRR_list.txt)
do
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/CIRI2/${i}/circ_candidates_convert.bed -o $working_dir/output/CIRI2/${i}/circularRNA_known.txt
done

# Merged
Rscript $script_dir/annotation.r -w $working_dir -a CIRI2 -i $working_dir/SRR_list.txt
