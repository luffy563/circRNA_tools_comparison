#!/usr/bin/env python
#-*- coding: UTF-8 -*-
import sys
import os
import pandas as pd
######################################################################################
syntax = '''
------------------------------------------------------------------------------------
Usage: python summarize.py segemehl_ouput_dir accession_id_file
Options:
    segemehl_ouput_dir: the directory of output for segemehl
    accession_id_file: the file contains accession_id of fastq files
Result: *.sum.bed summarize the raw results (_filt.sngl.bed) given by segemehl on all accession_id datasets 
------------------------------------------------------------------------------------
'''
######################################################################################

if len(sys.argv) != 3:
    print(syntax)
    sys.exit()
######################################################################################

## set the global variables
segemehl_ouput_dir = sys.argv[1]
accession_id_file = sys.argv[2]


f = open(accession_id_file, 'r')

for i in f.readlines():
    i = i.replace('\n', '')
    path = os.path.join(segemehl_ouput_dir,i,i+'_filt.sngl.bed')
    data = pd.read_csv(path,names=['chrom','start','end','name','score','strand'],sep='\t')
    data['counts'] = 1
    data = data.iloc[1:, :]
    data_1 = data.drop_duplicates(['chrom','start','end'])
    ID=[]
    for j in range(len(data_1['chrom'])):
        ID.append(j)
    print(len(ID))
    print(len(data_1['chrom']))
    data_1['ID'] = ID
    data_1.index = ID
    data_1 = data_1[['chrom','start','end','ID']]
    data_merged = pd.merge(data,data_1,on=['chrom','start','end'],how='left')
    data_merged = data_merged.drop(columns=['score'])
    data_2 = data_merged['counts'].groupby(data_merged['ID']).sum()
    data_1['counts']=data_2
    # data_sum = data.drop_duplicates()
    summarize_path = os.path.join(segemehl_ouput_dir,i,i+'.sum.bed')
    data_1.to_csv(summarize_path,index=False,sep = '\t',header=False)
    