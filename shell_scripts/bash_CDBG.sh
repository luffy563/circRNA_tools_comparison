##!/bin/bash

# CircDBG
for i in $(cat SRR_list.txt)
do
sed -ie 's/Reads1=/Reads1=SRA/${i}/${i}.fastq.gz/g' config.ini
./software/CircDBG-master/CircDBG/CircDBG/MakeFile/CircRNADBG ./config.ini
done

for i in $(cat SRR_list-1.txt)
do
sed -ie 's/Reads1=/Reads1=SRA/${i}/${i}_1.fastq.gz/g' config.ini
sed -ie 's/Reads2=/Reads2=SRA/${i}/${i}_2.fastq.gz/g' config.ini
./software/CircDBG-master/CircDBG/CircDBG/MakeFile/CircRNADBG ./config.ini
done


