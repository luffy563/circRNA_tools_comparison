#!/bin/bash

# path=software/KNIFE-master/circularRNApipeline_Standalone
# # create the exonDB
# path_to_output_directory=/home/data1/liuhongfei/pred_circ/ref_genome/KNIFE
# python makeExonDB.py -f genome.fasta -a genome.gtf -o path_to_output_directory
# # creates the junction indices
# path_to_circularRNApipeline=/home/data/liuhongfei/pred_circ/software/KNIFE-master/circularRNApipeline_Standalone
# sh createJunctionIndex.sh $path_to_circularRNApipeline $path_to_output_directory hg19 1000000

for i in $(cat SRR_list.txt)
do
mkdir -p output/KNIFE
read_directory=/home/data/liuhongfei/pred_circ/SRA/SRR9077/${i}
read_id_style=complete
alignment_parent_directory=output/KNIFE
dataset_name=${i}
junction_overlap=15
mode=sam
report_directory_name=circRNAreads
ntrim=50
denovoCircMode=1
junction_id_suffix=hsa_

cd $path
sh completeRun.sh read_directory read_id_style alignment_parent_directory dataset_name junction_overlap mode report_directory_name ntrim denovoCircMode junction_id_suffix 2>&1 | tee out.log
done

for i in $(cat SRR_list-1.txt)
do
mkdir -p output/KNIFE
read_directory=/home/data/liuhongfei/pred_circ/SRA/${i}/${i}
read_id_style=complete
alignment_parent_directory=output/KNIFE
dataset_name=${i}
junction_overlap=13
mode=sam
report_directory_name=circRNAreads
ntrim=66
denovoCircMode=1
junction_id_suffix=hsa_

cd $path
sh completeRun.sh read_directory read_id_style alignment_parent_directory dataset_name junction_overlap mode report_directory_name ntrim denovoCircMode junction_id_suffix 2>&1 | tee out.log
done