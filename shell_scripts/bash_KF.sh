#!/bin/bash

# KNIFE
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
script_dir=$working_dir/scripts
path=$working_dir/software/KNIFE-master/circularRNApipeline_Standalone

# create the exonDB
path_to_output_directory=$working_dir/KNIFE_exonDB
python $path/makeExonDB.py -f $working_dir/ref_genome/genome.fasta -a $working_dir/ref_genome/genome.gtf -o path_to_output_directory
# creates the junction indices
path_to_circularRNApipeline=$working_dir/software/KNIFE-master/circularRNApipeline_Standalone
sh $path/createJunctionIndex.sh $path_to_circularRNApipeline $path_to_output_directory hg19 1000000

# Alignment
simulated_datasets=(background pos mixed)
for i in $(cat $working_dir/SRR_list.txt)
do
mkdir -p $working_dir/output/KNIFE

if [[ "${simulated_datasets[@]}" =~ "$i" ]]
then echo $i
read_id_style=appended
else
read_id_style=complete
fi
read_directory=$working_dir/SRA/${i}/${i}
alignment_parent_directory=$working_dir/output/KNIFE
dataset_name=${i}
junction_overlap=15
mode=sam
report_directory_name=circRNAreads
ntrim=50
denovoCircMode=1
junction_id_suffix=hsa_

# cd $path
sh $path/completeRun.sh $read_directory $read_id_style $alignment_parent_directory $dataset_name $junction_overlap $mode $report_directory_name $ntrim $denovoCircMode $junction_id_suffix
done

## conversion 
# cd $working_dir
Rscript $script_dir/conversion.r -w $working_dir -a KF -i $working_dir/SRR_list.txt

# annotation
for i in $(cat $working_dir/SRR_list.txt)
do
# liftOver
liftOver $working_dir/output/KNIFE/${i}/circ_candidates_convert.bed $working_dir/ref_genome/liftover/hg19ToHg38.over.chain $working_dir/output/KNIFE/${i}/circ_candidates_convert_liftOver.bed $working_dir/output/KNIFE/${i}/unmapped.bed
# Annotation
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/KNIFE/${i}/circ_candidates_convert_liftOver.bed -o $working_dir/output/KNIFE/${i}/circularRNA_known.txt

done

Rscript $script_dir/annotation.r -w $working_dir -a KF -i $working_dir/SRR_list.txt
