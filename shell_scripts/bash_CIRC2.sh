#!/bin/bash

# CIRCexplorer2
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts

# Alignment (paired-end seq)
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p $working_dir/output/tophat/${i}
mkdir -p $working_dir/output/CIRCexplorer2/${i}
tophat2 -o $working_dir/output/CIRCexplorer2 -p 20 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $working_dir/bw_index/genome $working_dir/SRA/${i}/${i}_1.fastq $working_dir/SRA/${i}/${i}_2.fastq

# Parse
CIRCexplorer2 parse --pe -t TopHat-Fusion $working_dir/output/CIRCexplorer2/${i}/accepted_hits.bam > $working_dir/output/CIRCexplorer2/${i}/tophat_fusion.log
done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a CIRC2 -i $working_dir/SRR_list.txt

# Annotate
for i in $(cat $working_dir/SRR_list.txt)
do

CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/CIRCexplorer2/${i}/back_spliced_junction.bed -o $working_dir/output/CIRCexplorer2/${i}/circularRNA_known.txt
done

Rscript $script_dir/annotation.r -w $working_dir -a CIRC2 -i $working_dir/SRR_list.txt
