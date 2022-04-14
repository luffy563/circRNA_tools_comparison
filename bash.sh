#!/bin/bash
####################################################################
#  One-step script for running all software on different datasets  #
####################################################################
## Date: 2022-4-13
## Author: Hongfei Liu
## Contact: lhf563@nwafu.edu.cn

## The tree structure of working directory

## Generate background datasets
# ART version: 2.5.8
working_directory=$HOME # where you can change the working_directory
art=$working_directory/software/circRNA_detection_review-master/simu/art_illumina
## Change the simulated depths for background datasets (-f 100 -> -f 200 on 20210705)
for coverage in $(seq 5 5 30)
do
mkdir -p $working_directory/SRA/background/background_$coverage/
$art -ss HS25 -d background -na -i ./ref_mrna/ref_mrna.fa -o $working_directory/SRA/background_$coverage/background -l 101 -f $coverage -p -m 350 -s 10 -sp -rs 20210705 -qs -13 -qs2 -13
mv $working_directory/SRA/background_$coverage/background1.fq $working_directory/SRA/background_$coverage/background_${coverage}_1.fastq && mv $working_directory/SRA/background_$coverage/background2.fq $working_directory/SRA/background_$coverage/background_${coverage}_2.fastq && echo "done"
done

## Generate positive datasets
for coverage in $(seq 5 5 30)
do
mkdir -p $working_directory/SRA/postest/pos_$coverage/
perl $working_directory/software/circRNA_detection_review-master/simu/CIRIsimulator.pl -1 $working_directory/SRA/pos_$coverage/pos_${coverage}_1.fastq -2 $working_directory/SRA/pos_$coverage/pos_${coverage}_2.fastq -O $working_directory/SRA/pos_circRNAs.list -G $working_directory/ref_genome/hg38.refGene.gtf -DB $working_directory/ref_circ/combine_final_gh38_database.txt -C $coverage -LC 0 -R 1 -LR 1 -L 100 -E 10 -I 350 -D $working_directory/ref_genome/hg38 -CHR1 0 -M 50
done

## Combine positive and background datasets into mixed datasets
for coverage in $(seq 5 5 30)
do
mkdir -p $working_directory/SRA/mixed_$coverage/
cat $working_directory/SRA/background_$coverage/background_${coverage}_1.fastq $working_directory/SRA/pos_$coverage/pos_${coverage}_1.fastq > $working_directory/SRA/mixed_$coverage/mixed_${coverage}_1.fastq
cat $working_directory/SRA/background_$coverage/background_${coverage}_2.fastq $working_directory/SRA/pos_$coverage/pos_${coverage}_2.fastq > $working_directory/SRA/mixed_$coverage/mixed_${coverage}_2.fastq
done

## Get the reads length of all datasets
python get_reads_length.py $working_directory/SRA $working_directory/readsLen.txt


### Analhysis of software on datasets
## Create the log directory
mkdir -p $working_directory/output/logs/

running_dir=$working_directory/shell_scripts
# # find_circ
{ time nohup bash $running_dir/bash_FC.sh $working_directory &>output/logs/FC.log & } 2>output/logs/FC.time

# # CIRCexplorer2
{ time nohup bash $running_dir/bash_CIRC2.sh $working_directory &>output/logs/CIRC2.log & } 2>output/logs/CIRC2.time

# # CircRNAfinder
{ time nohup bash $running_dir/bash_CF.sh $working_directory &>output/logs/CF.log & } 2>output/logs/CF.time

# # CDBG
{ time nohup bash $running_dir/bash_CDBG.sh $working_directory &>output/logs/CDBG.log & } 2>output/logs/CDBG.time

# # CircMarker
{ time nohup bash $running_dir/bash_CM.sh $working_directory &>output/logs/CM.log & } 2>output/logs/CM.time

# # CIRI2
{ time nohup bash $running_dir/bash_CIRI2.sh $working_directory &>output/logs/CIRI2.log & } 2>output/logs/CIRI2.time

# # DCC
{ time nohup bash $running_dir/bash_DCC.sh $working_directory &>output/logs/DCC.log & } 2>output/logs/DCC.time

# # Mapsplice
{ time nohup bash $running_dir/bash_MP.sh $working_directory &>output/logs/MP.log & } 2>output/logs/MP.time

# # Segemeghl
{ time nohup bash $running_dir/bash_SE.sh $working_directory &>output/logs/SE.log & } 2>output/logs/SE.time

# # KNIFE
{ time nohup bash $running_dir/bash_KF.sh $working_directory &>output/logs/NCL.log & } 2>output/logs/NCL.time

# # UROBORUS
{ time nohup bash $running_dir/bash_UB.sh $working_directory &>output/logs/UB.log & } 2>output/logs/UB.time

