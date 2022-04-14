##!/bin/bash

# CircMarker
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
config_dir=$working_dir/config
script_dir=$working_dir/scripts
path=$working_dir/software/CircMarker-master/CircRnaDetectDraft/CircRNA_Detect_Draft/MakeFile

# Modify the config file and run analysis (paired-end seq)
cat $working_dir/readsLen.txt | while read i len
do
mkdir -p $working_dir/output/CircMarker/${i}
# cp $config_dir/CM_config.ini $path/CM_config.ini
# result=SRA\/${i}\/${i}.fastq
# echo $result
sed -i "/Reads1=/s/$/SRA\/${i}\/${i}_1.fastq/" $path/CM_config.ini
sed -i "/Reads2=/s/$/SRA\/${i}\/${i}_2.fastq/" $path/config.ini
sed -i "/ReadsLen=/s/$/$len/" $path/CM_config.ini
$path/CircRnaDetectDraft $path/CM_config.ini
cp -r $working_dir/Detection_Result $working_dir/output/CircMarker/${i}
done

# Conversion
Rscript $script_dir/conversion.r -w $working_dir -a CM -i $working_dir/SRR_list.txt

# Annotate
for i in $(cat temp.txt)
do
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/CircMarker/${i}/circ_candidates_convert.bed -o $working_dir/output/CircMarker/${i}/circularRNA_known.txt
done

# Annotation
Rscript $script_dir/annotation.r -w $working_dir -a CM -i $working_dir/SRR_list.txt


