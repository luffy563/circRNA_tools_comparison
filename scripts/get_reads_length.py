#!/usr/bin/env python
#-*- coding: UTF-8 -*-
import sys
import os
######################################################################################
syntax = '''
------------------------------------------------------------------------------------
Usage: python get_reads_lenght.py fastq_file_dir output_file
Options:
    fastq_file_dir: the directory of .fastq or .fq files
    output_file: the absolute or relative path of average read length for all fastq_file
Result: .txt file same name as input name plus "_lenghts.txt" 
------------------------------------------------------------------------------------
'''
######################################################################################

if len(sys.argv) != 3:
    print(syntax)
    sys.exit()

######################################################################################
## set the global variables
fastq_file_dir = sys.argv[1]
output_file = sys.argv[2]
## get the fastq files in fastq_file_dir
accession_ids = os.listdir(fastq_file_dir)

# outfile = open(prefix + '_' + 'lenghts.txt', 'w')
def get_average_readLen(fastq_file_path):
    fastq_file = open(fastq_file_path,'r')
    seq = ''
    readLenSum = 0
    readNumber = 0
    for line in fastq_file:
        line = line.rstrip('\n')
        if line.startswith('@'):
            if seq:
                readLenSum += len(seq)
                readNumber += 1
                seq = ""
            name = line 
        else:
            seq = line 
    
    fastq_file.close()
    return int(readLenSum/readNumber)


with open(output_file,'w') as f:
    for accession_id in accession_ids:
        fastq_files = os.listdir(os.path.join(fastq_file_dir,accession_id))
        print(fastq_files)
        prefix = fastq_files[0].split('.')[0]
        suffix = fastq_files[0].split('.')[1]
        if suffix not in ['fastq','fq']:
            print("Compressed format is not supported!")
            break
        
        fastq_file = prefix + "." + suffix
        fastq_file_path = os.path.join(fastq_file_dir,accession_id,fastq_file)
        
        readAveLen = get_average_readLen(fastq_file_path)
        print(readAveLen)
        f.write(accession_id+"\t"+str(readAveLen)+"\n")