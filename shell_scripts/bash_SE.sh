#!/bin/bash

# segemehl
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts

mkdir -p $working_dir/se_index/
# Indexing
segemehl.x -x $working_dir/se_index/gh38.idx -d $working_dir/ref_genome/genome.fa
# single-end seq
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p $working_dir/output/segemehl/${i}
# Alignment paired-end seq
segemehl.x -t 20 -i $working_dir/se_index/gh38.idx -d $working_dir/ref_genome/genome.fa -q $working_dir/SRA/${i}/${i}_1.fastq $working_dir/SRA/${i}/${i}_2.fastq -S -o $working_dir/output/segemehl/${i}
done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a SE -i $working_dir/SRR_list.txt
for i in $(cat $working_dir/SRR_list.txt)
do
# Annotate
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/segemehl/${i}/circ_candidates-1.bed -o $working_dir/output/segemehl/${i}/circularRNA_known.txt
done
Rscript $script_dir/annotation.r -w $working_dir -a SE -i $working_dir/SRR_list.txt

