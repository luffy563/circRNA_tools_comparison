#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


f = open('software_list.txt', 'r') # read software tools list
background_file = open('background.txt', 'r') # read background dataset id

softlist, backgrounds = [], []

for soft in f:
    softlist.append(soft.replace('\n',''))

for background in background_file:
    backgrounds.append(background.replace('\n',''))

# print(softlist)
# print(backgrounds)

colors = ['#CEBA89',
 '#8D483E',
 '#F2FBAE',
 '#18B3AA',
 '#9F9A48',
 '#FB989B',
 '#55C98E',
 '#1A5CC8',
 '#AABD57',
 '#46EFFD',
 '#6C473C']

background_datasets = {}
for back in backgrounds:
    background_datasets[back] = []
    for line in softlist:
        # line = line.replace('\n', '')
        path = './' + line + '/' + back + '/circ_candidates_convert.bed'
        if line == 'find_circ':
            # path = './' + line + '/' + pos + '/circ_candidates-1.bed'
            data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'n_reads',\
                                              'strand', 'n_uniq', 'uniq_bridges', 'best_qual_left',\
                                              'best_qual_right', 'tissues', 'tiss_counts', 'edits',\
                                              'anchor_overlap', 'breakpoints', 'signal', 'strandmatch', 'category'])
            # data['chrom'] = pd.Series([('chr'+x) for x in data['chrom']])
            data_name = data['name'].str.split('/', expand=True).iloc[:, [0, 1]]
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        elif line == 'CIRCexplorer2':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand','unknown'])
            data_name = data['name'].str.split('/', expand=True).iloc[:,[0,1]]
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)

        elif line == 'CircRNAfinder':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)

        elif line == 'Mapsplice':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path).iloc[:,0:4]
            data.columns = ['chrom', 'start', 'end', 'name']
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        elif line == 'CIRI2':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            # if (data)
            try:
                data = pd.read_table(path).iloc[:, 0:6]
                data.columns = ['chrom', 'start', 'end', 'name', 'score', 'linear']
                data_name = data['name'].str.split('/', expand=True)
                data[['name', 'n_reads']] = data_name
                data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
                background_datasets[back].append(data)

            except:
                data = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'linear'])
                data.columns = ['chrom', 'start', 'end', 'name', 'score', 'linear']
                # data_name = data['name'].str.split('/', expand=True)
                # data[['name', 'n_reads']] =
                data['n_reads'] = 'NA'
                background_datasets[back].append(data)

        elif line == 'NCLscan':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path)
            data.columns = ['chrom', 'start',  'end', 'name']
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        elif line == 'DCC':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path).iloc[:, 0:4]
            data.columns = ['chrom', 'start',  'end', 'name']
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        elif line == 'UROBORUS':
            path = './' + line + '/' + back + '/remapped_circ_candidates_convert_remapped.bed'
            print(path)
            data = pd.read_table(path)
            data.columns = ['chrom', 'start',  'end', 'name']
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        elif line == 'segemehl':
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path)
            data.columns = ['chrom', 'start', 'end', 'name', 'score']
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        elif line == 'KNIFE':
            path = './' + line + '/' + back + '/circ_candidates_convert_liftOver.bed'
            data = pd.read_table(path,header=None)
            data.columns = ['chrom', 'start', 'end', 'name', 'score']
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        else:
            # path = './' + line + '/pos/circ_candidates_convert.bed'
            data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'strand', 'type', 'score'])
            data_name = data['name'].str.split('/', expand=True)
            data[['name', 'n_reads']] = data_name
            data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
            background_datasets[back].append(data)
        # path = './' + soft + '/' + back + '/circ_candidates_convert.bed'

# print(background_datasets)
# print(len(background_datasets))
# print(len(background_datasets['background']))

## Fig - Upset plot of background datasets for different software
cols = ['datasets','software','circRNA_id','readcounts']
upset_df = {}.fromkeys(cols,[])
for back in backgrounds:
    # each_soft_data = background_datasets[back][0]
    # upset_df['backgrounds'] += [back] * len(each_soft_data['chrom'])
    for soft in softlist:
        each_soft_data = background_datasets[back][softlist.index(soft)]
        upset_df['datasets'] = upset_df['datasets'] + [back] * len(each_soft_data['chrom'])
        upset_df['software'] = upset_df['software'] + [soft] * len(each_soft_data['chrom'])
        upset_df['circRNA_id'] = upset_df['circRNA_id'] + [str(chrom)+":"+str(start)+"-"+str(end)
                                                           for chrom, start, end in zip(each_soft_data['chrom'],each_soft_data['start'],
                                                                       each_soft_data['end'])]
        upset_df['readcounts'] = upset_df['readcounts'] + [i for i in each_soft_data['n_reads']]

## Output the raw expression dataset of sofware on background datasets
upset_df = pd.DataFrame(upset_df)
# print(upset_df)
upset_df.to_csv('./BackgroundRaw.csv',sep=',',index=False,header=True)

## Output the dataset for plotting the Upset Plot
upset_df.drop('readcounts',inplace=True,axis=1)
upset_df_table = upset_df.pivot_table(index='circRNA_id',columns='software',aggfunc='count',fill_value=0)
upset_df_table.columns = [col[1] for col in upset_df_table.columns.tolist()]
upset_df_table.to_csv('./BackgroundPlotdf.csv',header=True,index=True)
