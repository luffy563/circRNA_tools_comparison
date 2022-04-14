#!/bin/bash

# Mapsplice
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts
path=$working_dir/Mapsplice/bin

# Alignment
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p output/Mapsplice/$i
python $path/mapsplice.py -p 20 -c $working_dir/ref_genome -x $working_dir/bw_index/genome -1 $working_dir/SRA/$i/${i}_1.fastq -o $working_dir/output/Mapsplice/${i}_1 --min-fusion-distance 200 --gene-gtf $working_dir/ref_genome/gh38.gtf
done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a MP -i $working_dir/SRR_list.txt

# Annotate
for i in $(cat $working_dir/SRR_list.txt)
do
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/Mapsplice/${i}/circ_candidates_convert.bed -o $working_dir/output/Mapsplice/${i}/circularRNA_known.txt
done

# Merged
Rscript $script_dir/annotation.r -w $working_dir -a MP -i $working_dir/SRR_list.txt