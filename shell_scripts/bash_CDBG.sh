##!/bin/bash

# CircDBG
echo "Working directory:" $1;
# Set the global varibales
working_dir=$1
config_dir=$working_dir/config
script_dir=$working_dir/scripts
# Modify the config file and run analysis
cat $working_dir/readsLen.txt | while read i len
do
mkdir -p $working_dir/output/CircDBG/

sed -i 's/Reads1=/s/$/SRA\/${i}\/${i}_1.fastq/' $config_dir/CDBG_config.ini
sed -i 's/Reads2=/s/$/SRA/\${i}\/${i}_2.fastq/' $config_dir/CDBG_config.ini
sed -i 's/ReadLen=/s/$/$len/g' $config_dir/CDBG_config.ini
$working_dir/software/CircDBG-master/CircDBG/CircDBG/MakeFile/CircRNADBG $config_dir/CDBG_config.ini
done

# Conversion raw format to bed format (circ_candidates_convert.bed)
Rscript $script_dir/conversion.r -w $working_dir -a CDBG -i $working_dir/SRR_list.txt

# Annotation by CIRCexplorer2 annotate moudle
for i in $(cat $working_dir/SRR_list.txt)
do
CIRCexplorer2 annotate -r $working_dir/ref_genome/gh38_ref.txt -g $working_dir/ref_genome/gh38_chrID.fa -b $working_dir/output/CircDBG/${i}/circ_candidates_convert.bed -o $working_dir/output/CircDBG/${i}/circularRNA_known.txt
done

# Merged by combined_final_hg38.txt
Rscript $script_dir/annotation.r -w $working_dir -a CDBG -i $working_dir/SRR_list.txt
