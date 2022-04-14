#!/bin/bash

# find_circ
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts

# Alignment for paired-end seq
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir $working_dir/output/find_circ/${i}
# Alignment
bowtie2 -p 8 --very-sensitive --score-min=C,-15,0 --mm -x $working_dir/bw2_index/genome -q -1 $working_dir/SRA/${i}/${i}_1.fastq -2 $working_dir/SRA/${i}/${i}_2.fastq | samtools view -hbuS - | samtools sort -o $working_dir/output/find_circ/${i}/output.bam
# Filter
samtools view -hf 4 $working_dir/output/find_circ/${i}/output.bam | samtools view -Sb - > $working_dir/output/find_circ/${i}/unmapped.bam
python $working_dir/software/find_circ/find_circ-1.2/unmapped2anchors.py $working_dir/output/find_circ/${i}/unmapped.bam | gzip > $working_dir/output/find_circ/${i}/anchors.fastq.gz
# Realignment
bowtie2 -p 8 --reorder --mm --score-min=C,-15,0 -q -x $working_dir/bw2_index/genome -U $working_dir/output/find_circ/${i}/anchors.fastq.gz | python $working_dir/software/find_circ/find_circ-1.2/find_circ.py -G $working_dir/ref_genome/genome.fa -p hsa_ -s $working_dir/output/find_circ/${i}/find_circ.sites.log > $working_dir/output/find_circ/${i}/find_circ.sites.bed 2> $working_dir/output/find_circ/${i}/find_circ.sites.reads

# Filteration
grep CIRCULAR $working_dir/output/find_circ/${i}/find_circ.sites.bed | grep -v chrM | awk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE > $working_dir/output/find_circ/${i}/circ_candidates.bed
done

# conversion
Rscript $script_dir/conversion.r -w $working_dir -a FC -i $working_dir/SRR_list.txt
for i in $(cat $working_dir/SRR_list.txt)
do
# Annotate
CIRCexplorer2 annotate -r $working_dirref_genome/gh38_ref.txt -g $working_dirref_genome/gh38_chrID.fa -b $working_diroutput/find_circ/${i}/circ_candidates-1.bed -o $working_diroutput/find_circ/${i}/circularRNA_known.txt
done
Rscript $script_dir/annotation.r -w $working_dir -a FC -i $working_dir/SRR_list.txt

