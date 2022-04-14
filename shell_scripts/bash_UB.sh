#!/bin/bash

# UROBORUS 
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts

for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p $working_dir/output/tophat/${i}
mkdir -p $working_dir/output/UROBORUS/$i
rm -rf $working_dir/circRNA_list.txt
tophat -p 20 -o $working_dir/output/tophat/${i} $working_dir/bw2_index/genome $working_dir/SRA/${i}/${i}_1.fastq $working_dir/SRA/${i}/${i}_2.fastq
samtools view $working_dir/output/tophat/${i}/unmapped.bam > $working_dir/output/tophat/${i}/unmapped.sam
perl $working_dir/software/UROBORUS-master/bin/UROBORUS.pl -p 20 -index $working_dir/bw_index/hg19/genome \
-gtf $working_dir/ref_genome/genes.gtf -fasta $working_dir/ref_genome/hg19 \
$working_dir/output/tophat/$i/unmapped.sam output/tophat/$i/accepted_hits.bam
mv $working_dir/circRNA_list.txt $working_dir/output/UROBORUS/$i/circRNA_list.txt

done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a UB -i $working_dir/SRR_list.txt

# remap by liftover
for i in $(cat $working_dir/SRR_list.txt)
do
liftOver $working_dir/output/UROBORUS/$i/circ_candidates_convert_remapped.bed $working_dir/remap/hg19ToHg38.over.chain \
$working_dir/output/UROBORUS/$i/remapped_circ_candidates_convert_remapped.bed $working_dir/output/UROBORUS/$i/remapped_circ_candidates_convert_remapped.bed

done

# Annotate
for i in $(cat $working_dir/SRR_list.txt)
do
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/UROBORUS/${i}/remapped_circ_candidates_convert_remapped.bed -o $working_dir/output/UROBORUS/${i}/circularRNA_known.txt
done

# Merged
Rscript $script_dir/annotation.r -w $working_dir -a UB -i $working_dir/SRR_list.txt
