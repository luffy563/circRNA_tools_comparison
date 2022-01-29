##!/bin/bash

# CircMarker

path=software/CircMarker-master/CircRnaDetectDraft/CircRNA_Detect_Draft/MakeFile
for i in $(cat SRR_list.txt)
do
mkdir -p output/CircMarker/${i}
cp config-1.ini $path/config.ini
result=SRA\/SRR9077\/${i}.fastq.gz
echo $result
sed -i "/Reads1=/s/$/SRA\/SRR9077\/${i}.fastq.gz/" $path/config.ini
sed -i "/ReadsLen=/s/$/76/" $path/config.ini
./$path/CircRnaDetectDraft ./$path/config.ini
cp -r Detection_Result output/CircMarker/${i}
done


# paired-end seq
for i in $(cat SRR_list-1.txt)
do
mkdir -p output/CircMarker/${i}
cp config-1.ini $path/config.ini
result=SRA\/${i}\/${i}_1/2.fastq.gz
echo $result
sed -i "/Reads1=/s/$/SRA\/${i}\/${i}_1.fastq.gz/" $path/config.ini
sed -i "/Reads2=/s/$/SRA\/${i}\/${i}_2.fastq.gz/" $path/config.ini
sed -i "/ReadsLen=/s/$/100/" $path/config.ini
./$path/CircRnaDetectDraft ./$path/config.ini
cp -r Detection_Result output/CircMarker/${i}
done


