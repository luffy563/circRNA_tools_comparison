#!/bin/bash

#Mapsplice
path=~/miniconda3/envs/py2/bin

for i in $(cat SRR_list.txt)
do
mkdir -p output/Mapsplice/$i
python $path/mapsplice.py -p 16 -c ref_genome -x bw_index/genome -1 SRA/$i/$i.fastq.gz -o output/Mapsplice/$i --min-fusion-distance 200 --gene-gtf ref_genome/gh38.gtf
done