import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from sklearn.metrics import auc

# # function for selecting color set
# def randomcolor(n):
    # colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    # color = ""
    # colors=[]
    # for i in range(n):
        # color=''
        # for j in range(6):
            # color += colorArr[random.randint(0,14)]
        # color = '#' + color
        # colors.append(color)
    # return colors

# colors = randomcolor(11)
# plt.bar(range(0,9),1,color=colors)
# plt.show()

# Set the fixed color set for different software
colors = ['#f8766d','#d39200','#93aa00','#00ba38','#00c19f','#00b9e3','#619cff','#db72fb','#ff61c3']


# total reads counts
path = './' + 'CIRI2' + '/positive/circ_candidates_convert.bed'
data = pd.read_table(path).iloc[:, 0:6]
data.columns = ['chrom', 'start', 'end', 'name', 'score', 'linear']
data_name = data['name'].str.split('/', expand=True)
data[['name', 'n_reads']] = data_name
data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
data['total_counts'] = data['n_reads'] + data['linear']
data_total_counts = data[['chrom', 'start', 'end', 'total_counts']]
data_total_counts.to_csv('data_total_counts.csv',index=False, sep=',')
## Positive dataset
############
# Figure 1 #
############
recovery_matrix = pd.DataFrame(data=None)

fig=plt.figure(figsize=(12,8),dpi = 300)
P_pos, R_pos, AUC_pos, softlist = [], [], [], []
plt.style.use('ggplot')
for line,color in zip(f, colors):
    line = line.replace('\n', '')
    softlist.append(line)
    if line == 'find_circ':
        # path = './' + line + '/positive/circ_candidates-1.bed'
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'n_reads',\
                                          'strand', 'n_uniq', 'uniq_bridges', 'best_qual_left',\
                                          'best_qual_right', 'tissues', 'tiss_counts', 'edits',\
                                          'anchor_overlap', 'breakpoints', 'signal', 'strandmatch', 'category'])
        # data['chrom'] = pd.Series([('chr'+x) for x in data['chrom']])
    elif line == 'CIRCexplorer2':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand','unknown'])
        data_name = data['name'].str.split('/', expand=True).iloc[:,[0,1]]
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])


    elif line == 'CircRNAfinder':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])

    elif line == 'Mapsplice':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path).iloc[:, 0:6]
        data.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'CIRI2':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path).iloc[:, 0:6]
        data.columns = ['chrom', 'start', 'end', 'name', 'score', 'linear']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])

    elif line == 'NCLscan':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path)
        data.columns = ['chrom', 'start',  'end', 'name']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'DCC':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path).iloc[:, 0:4]
        data.columns = ['chrom', 'start',  'end', 'name']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'UROBORUS':
        path = './' + line + '/pos/remapped_circ_candidates_convert_remapped.bed'
        data = pd.read_table(path)
        data.columns = ['chrom', 'start',  'end', 'name']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'segemehl':
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path)
        data.columns = ['chrom', 'start', 'end', 'name', 'score']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])

    else:
        path = './' + line + '/pos/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'strand', 'type', 'score'])
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])


    tp_path = './' + line + '/positive/circRNA_known_annotated.csv'
    print(tp_path)
    tp_data = pd.read_csv(tp_path, header=0)
    # recovered_data = tp_data[['exonSizes','readNumber']]
    #
    # recovered_data_exonsize = recovered_data['exonSizes'].str.split(',', expand=True)
    # recovered_data_exonsize = recovered_data_exonsize.fillna(value=0)
    # recovered_data_exonsize=recovered_data_exonsize.astype('int')
    # recovered_data_exonsize = np.sum(recovered_data_exonsize,axis=1)
    # recovered_data['length'] = recovered_data_exonsize
    # recovered_data = recovered_data[(recovered_data['readNumber']>0) & (recovered_data['length']>350)]
    # recovery_rate = recovered_data['readNumber']*100/recovered_data['length']
    #
    # recovery_matrix[line] = recovery_rate/(recovery_rate.max()-recovery_rate.min())

    # recovery rate
    recovered_data = pd.merge(tp_data, data_total_counts,on=['chrom','start','end'])
    recovery_matrix[line] = recovered_data['readNumber']/recovered_data['total_counts']
    data = data[data['n_reads'] >=2]
    precision, recall = [],[]
    pd_data = tp_data[['chrom','start','end','circRNA_ID']]
    cross_data = pd.merge(data,pd_data,how='left',on=['chrom','start','end'])
    cross_data = cross_data.drop_duplicates(subset=['chrom','start','end'], keep='first')
    # for i in range(len(recovered_data)):
    # cross_data = cross_data[cross_data['n_reads'] >= 10]
    for i in range(max(data['n_reads']), min(data['n_reads'])-1, -1):
        TP = len(cross_data[(cross_data['n_reads'] >= i) & (pd.notnull(cross_data['circRNA_ID']))])
        P = TP/(len(cross_data[cross_data['n_reads'] >= i])+1)
        R = TP/(len(cross_data[pd.notnull(cross_data['circRNA_ID'])])+1)
        precision.append(P)
        recall.append(R)
    precision = np.array(precision)
    recall = np.array(recall)
    auc_value_pos = auc(recall,precision)
    AUC_pos.append(auc_value_pos)
    P_pos.append(precision[-2])
    R_pos.append(recall[-2])

    plt.plot(recall, precision, c=color, label=line,linewidth=4.0)
    plt.scatter(recall, precision, c=color, s=100)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.grid(b=False)
    plt.tick_params(labelsize=18)
    plt.legend(loc='lower right',fontsize=15)
    plt.xlabel('Recall',fontsize=25)
    plt.ylabel('Precision',fontsize=25)

plt.savefig('fig1.svg',dpi=600)
plt.show()


#############
# Figure S4 #
#############
# figure S4 recoverd rate at each candidate for different software
fig=plt.figure(figsize=(12,8),dpi = 300)
plt.style.use('ggplot')
# recovery_matrix_1=pd.melt(recovery_matrix,var_name='Software',value_name='Recovery rate')
ax = recovery_matrix.boxplot(fontsize=20,patch_artist=True,boxprops={'color': 'orangered', 'facecolor': 'blue'}, showfliers=False,medianprops={'linestyle':'-','color':'black'})
ax.grid(b=False)
recovery_matrix.to_csv('recovery_matrix.csv',index=False,sep=',')
for tick in ax.get_xticklabels():
    tick.set_rotation(45)
plt.savefig('figS4.svg',dpi=600)
plt.show()


### Mixed datasets
############
# Figure 2 #
############
# Figure 2
f = open('software_list.txt', 'r') # read software list
fig=plt.figure(figsize=(12, 8), dpi=300)
plt.style.use('ggplot')
subset = []
P_mix, R_mix, AUC_mix = [], [], []
for line,color in zip(f, colors):
    line = line.replace('\n', '')
    if line == 'find_circ':
        path = './' + line + '/Mixed/circ_candidates-1.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'n_reads',\
                                          'strand', 'n_uniq', 'uniq_bridges', 'best_qual_left',\
                                          'best_qual_right', 'tissues', 'tiss_counts', 'edits',\
                                          'anchor_overlap', 'breakpoints', 'signal', 'strandmatch', 'category'])
        # data['chrom'] = pd.Series([('chr'+x) for x in data['chrom']])

    elif line == 'CIRCexplorer2':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand','unknown'])
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])

        #data[['chrom', 'start', 'end', 'name', 'n_reads','score','strand','unknown']] = data
   # elif line == 'CircDBG':
   #      path = './' + line + '/positive/Detection_Result/back_spliced_junction-1.bed'
   #      data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score','strand',])
    elif line == 'CircRNAfinder':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'Mapsplice':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path).iloc[:, 0:6]
        data.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'CIRI2':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path).iloc[:, 0:5]
        data.columns = ['chrom', 'start', 'end', 'name', 'score']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    # print(path)
    elif line == 'NCLscan':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path)
        data.columns = ['chrom', 'start',  'end', 'name']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'DCC':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path).iloc[:,0:4]
        data.columns = ['chrom', 'start',  'end', 'name']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'UROBORUS':
        path = './' + line + '/Mixed/remapped_circ_candidates_convert_remapped.bed'
        data = pd.read_table(path)
        data.columns = ['chrom', 'start',  'end', 'name']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    elif line == 'segemehl':
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path)
        data.columns = ['chrom', 'start', 'end', 'name', 'score']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])
    else:
        path = './' + line + '/Mixed/circ_candidates_convert.bed'
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'strand', 'type', 'score'])
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data['n_reads'] = pd.Series([int(x) for x in data['n_reads']])

    tp_path = './' + line + '/Mixed/circRNA_known_annotated.csv'
    print(tp_path)
    tp_data = pd.read_csv(tp_path, header=0)

    subset.append(tp_data['circRNA_ID'].tolist())
    precision, recall = [], []
    pd_data = tp_data[['chrom','start','end','circRNA_ID']]
    cross_data = pd.merge(data,pd_data,how='left',on=['chrom','start','end'])
    cross_data = cross_data.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
    cross_data = cross_data[cross_data['n_reads'] >= 2]

    for i in range(max(data['n_reads']), min(data['n_reads'])-1, -1):
        TP = len(cross_data[(cross_data['n_reads'] >= i) & (pd.notnull(cross_data['circRNA_ID']))])
        P = TP/len(cross_data[cross_data['n_reads'] >= i])
        R = TP/(len(cross_data[pd.notnull(cross_data['circRNA_ID'])])+1)
        precision.append(P)
        recall.append(R)
    precision = np.array(precision)
    recall = np.array(recall)
    auc_value_mix = auc(recall,precision)
    AUC_mix.append(auc_value_mix)
    P_mix.append(precision[-2]), R_mix.append(recall[-2])
    plt.plot(recall, precision, c=color, label=line, linewidth=4.0)
    plt.scatter(recall, precision, c=color, s=100)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.grid(b=False)
    plt.tick_params(labelsize=18)
    plt.legend(fontsize=15)
    plt.xlabel('Recall', fontsize=25)
    plt.ylabel('Precision', fontsize=25)
plt.savefig('fig1B.svg',dpi=600)
plt.show()

##############
# Figure 2-2 #
##############
# figure 2-2 comparison of PR value on different datasets
# fig2_2 = plt.figure(figsize=(12,8),dpi=300)

# plt.scatter(R_pos, P_pos, marker='o', c=colors, s=200)
# plt.scatter(R_mix, P_mix, marker='^', c=colors, s=200)
# ax = fig2_2.gca()
# handles,labels = ax.get_legend_handles_labels()
# plt.legend(handles, labels = softlist, loc='upper right')
#
# plt.show()

# The change of PR value between mixed and positive

### Real dataset
import matplotlib.pyplot as plt
import seaborn as sns

typelist=['reads_map','reads_map','reads_map','k-mer','k-mer','reads_map','reads_map','reads_map','reads_map']
hela_norm,hela_rnaser,k562_norm,k562_rnaser,yLabel = [],[],[],[],[]

# tp_data
# tp_path = './remapped_hsa_hg19_hela_circRNA.csv'
# tp_data = pd.read_csv(tp_path)
tp_path = 'remapped_hsa_hg19_circRNA.bed'
tp_data = pd.read_table(tp_path,header=None).iloc[:,0:4]
tp_data.columns = ['chrom','start','end','circRNA_ID']


# fig,ax=plt.subplots(nrows=2,ncols=1,figsize=(12,16),dpi=300)
# ax=ax.flatten()
# plt.style.use('ggplot')

for line,color in zip(softlist, colors):
    line = line.replace('\n', '')
    dataset,dataset_score=[],[]
    data_list = open('SRR_list-1.txt')
    precision, recall = [], []
    k562_precision, k562_recall = [], []
    for i in data_list:
        i = i.replace('\n', '')
        path = './' + line + '/' + i + '/circularRNA_known.txt'
        data = pd.read_table(path, names=['chrom','start','end','name','score','strand','thickStart',\
            'thickEnd','itemRgb','exonCount','exonSizes','exonOffsets','readNumber',\
            'circType','geneName','isoformName','index','flankIntron'])
        # data['chrom'] = pd.Series([('chr'+x) for x in data['chrom']])
        # data = data[['chrom','start','end','name','readNumber','isoformName']]
        data = data[['chrom', 'start', 'end', 'readNumber']]
        # data_score = data[['score']]
        # dataset_score.append(data_score)
        dataset.append(data)

    hela_norm_data=pd.merge(dataset[0],dataset[1],how='outer',on=['chrom','start','end'])
    hela_rnaser_data=pd.merge(dataset[2],dataset[3],how='outer',on=['chrom','start','end'])
    k562_norm_data=pd.merge(dataset[4],dataset[5],how='outer',on=['chrom','start','end'])
    k562_rnaser_data=pd.merge(dataset[6],dataset[7],how='outer',on=['chrom','start','end'])

    # Mean value of before and after RNaseR treat
    norm_data=pd.merge(hela_norm_data,hela_rnaser_data,how='outer',on=['chrom','start','end'])
    hypo_data=pd.merge(k562_norm_data,k562_rnaser_data,how='outer',on=['chrom','start','end'])
    # PR-curve
    hela_norm_data = hela_norm_data.fillna(0)
    hela_rnaser_data = hela_rnaser_data.fillna(0)
    k562_norm_data = k562_norm_data.fillna(0)
    k562_rnaser_data = k562_rnaser_data.fillna(0)

    norm_data=norm_data.fillna(0)
    hypo_data = hypo_data.fillna(0)

    norm_data['read_counts'] = (norm_data['readNumber_x_x'] + norm_data['readNumber_x_y'] + norm_data['readNumber_y_x'] + norm_data['readNumber_y_y'])/4
    hypo_data['read_counts'] = (hypo_data['readNumber_x_x'] + hypo_data['readNumber_x_y'] + hypo_data['readNumber_y_x'] + hypo_data['readNumber_y_y'])/4


    norm_cross_data = pd.merge(norm_data, tp_data, how='left', on=['chrom', 'start', 'end'])
    norm_cross_data = norm_cross_data.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
    hypo_cross_data = pd.merge(hypo_data, tp_data, how='left', on=['chrom', 'start', 'end'])
    hypo_cross_data = hypo_cross_data.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
    # for i in range(len(recovered_data)):\


    norm_cross_data = norm_cross_data[norm_cross_data['read_counts'] >= 2]
    hypo_cross_data = hypo_cross_data[hypo_cross_data['read_counts'] >= 2]
    # hela_dataset
    for i in range(round(max(norm_data['read_counts'])), round(min(norm_data['read_counts'])) - 1, -1):
        TP = len(norm_cross_data[(norm_cross_data['read_counts'] >= i) & (pd.notnull(norm_cross_data['circRNA_ID']))])
        P = TP / (len(norm_cross_data[norm_cross_data['read_counts'] >= i]) + 1)
        R = TP / (len(norm_cross_data[pd.notnull(norm_cross_data['circRNA_ID'])])+1)
        precision.append(P)
        recall.append(R)
    precision = np.array(precision)
    recall = np.array(recall)
    auc_value_pos = auc(recall, precision)
    AUC_pos.append(auc_value_pos)
    P_pos.append(precision[-2])
    R_pos.append(recall[-2])

    ax[0].plot(recall, precision, c=color, label=line, linewidth=4.0)
    ax[0].scatter(recall, precision, c=color, s=100)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    ax[0].grid(b=False)
    ax[0].tick_params(labelsize=15)
    ax[0].legend(fontsize=15)
    ax[0].set_xlabel('Recall', fontsize=20)
    ax[0].set_ylabel('Precision', fontsize=20)

    # k562_dataset
    for i in range(round(max(hypo_data['read_counts'])), round(min(hypo_data['read_counts'])) - 1, -1):
        TP = len(hypo_cross_data[(hypo_cross_data['read_counts'] >= i) & (pd.notnull(hypo_cross_data['circRNA_ID']))])
        P = TP / (len(hypo_cross_data[hypo_cross_data['read_counts'] >= i])+1)
        R = TP / (len(hypo_cross_data[pd.notnull(hypo_cross_data['circRNA_ID'])])+1)
        k562_precision.append(P)
        k562_recall.append(R)
    k562_precision = np.array(precision)
    k562_recall = np.array(recall)
    # auc_value_pos = auc(recall, precision)
    # AUC_pos.append(auc_value_pos)
    # P_pos.append(precision[-2])
    # R_pos.append(recall[-2])

    ax[1].plot(k562_recall, k562_precision, c=color, label=line, linewidth=4.0)
    ax[1].scatter(k562_recall, k562_precision, c=color, s=100)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.grid(b=False)
    plt.tick_params(labelsize=15)
    plt.legend(fontsize=15)
    plt.xlabel('Recall', fontsize=20)
    plt.ylabel('Precision', fontsize=20)

    hela_norm_data_sum=hela_norm_data[['readNumber_y', 'readNumber_x']]
    hela_rnaser_data_sum=hela_rnaser_data[['readNumber_y', 'readNumber_x']]
    k562_norm_data_sum = k562_norm_data[['readNumber_y', 'readNumber_x']]
    k562_rnaser_data_sum = k562_rnaser_data[['readNumber_y', 'readNumber_x']]

    hela_norm_data['n_reads']=hela_norm_data_sum.mean(axis=1)*1000000/sum(hela_norm_data_sum.mean(axis=1))
    hela_rnaser_data['n_reads']=hela_rnaser_data_sum.mean(axis=1)*1000000/sum(hela_rnaser_data_sum.mean(axis=1))
    k562_norm_data['n_reads'] = k562_norm_data_sum.mean(axis=1) * 1000000 / sum(k562_norm_data_sum.mean(axis=1))
    k562_rnaser_data['n_reads'] = k562_rnaser_data_sum.mean(axis=1) * 1000000 / sum(k562_rnaser_data_sum.mean(axis=1))

    hela_norm_data = hela_norm_data[['chrom', 'start', 'end', 'n_reads']]
    hela_rnaser_data = hela_rnaser_data[['chrom', 'start', 'end', 'n_reads']]
    k562_norm_data = k562_norm_data[['chrom', 'start', 'end', 'n_reads']]
    k562_rnaser_data = k562_rnaser_data[['chrom', 'start', 'end', 'n_reads']]

    hela_norm.append(hela_norm_data)
    hela_rnaser.append(hela_rnaser_data)
    k562_norm.append(k562_norm_data)
    k562_rnaser.append(k562_rnaser_data)
    yLabel.append(line)

# plt.savefig('fig1CD.svg',dpi=600)
# plt.show()

# norm_inter = pd.concat(norm, axis=1, join='inner')
# hypo_inter = pd.concat(hypo, axis=1, join='inner')
hela_norm_inter,hela_rnaser_inter,k562_norm_inter,k562_rnaser_inter = hela_norm[0],hela_rnaser[0],k562_norm[0],k562_rnaser[0]
for i in range(len(yLabel)-1):
    hela_norm_inter = pd.merge(hela_norm_inter,hela_norm[i+1],on=['chrom', 'start', 'end'])
    hela_rnaser_inter = pd.merge(hela_rnaser_inter, hela_rnaser[i + 1], on=['chrom', 'start', 'end'])
    k562_norm_inter = pd.merge(k562_norm_inter, k562_norm[i + 1], on=['chrom', 'start', 'end'])
    k562_rnaser_inter = pd.merge(k562_rnaser_inter, k562_rnaser[i + 1], on=['chrom', 'start', 'end'])
hela_norm_map = np.zeros((len(hela_norm_inter),len(yLabel)))
hela_rnaser_map = np.zeros((len(hela_rnaser_inter),len(yLabel)))
k562_norm_map = np.zeros((len(k562_norm_inter),len(yLabel)))
k562_rnaser_map = np.zeros((len(k562_rnaser_inter),len(yLabel)))

for i in range(len(yLabel)):
    hela_norm_map[:, i] = np.log2(hela_norm_inter.iloc[:, i + 3])
    hela_rnaser_map[:, i] = np.log2(hela_rnaser_inter.iloc[:, i + 3])
    k562_norm_map[:, i] = np.log2(k562_norm_inter.iloc[:, i + 3])
    k562_rnaser_map[:, i] = np.log2(k562_rnaser_inter.iloc[:, i + 3])
# Ordered by the software list
# hela_norm_map = hela_norm_map[np.lexsort(-hela_norm_map.T)].T
# hela_rnaser_map = hela_rnaser_map[np.lexsort(-hela_rnaser_map.T)].T
# k562_norm_map = k562_norm_map[np.lexsort(-k562_norm_map.T)].T
# k562_rnaser_map = k562_rnaser_map[np.lexsort(-k562_rnaser_map.T)].T
# Ranked by the mean value
hela_norm_map = hela_norm_map[np.argsort(-np.mean(hela_norm_map.T,axis=0))].T
hela_rnaser_map = hela_rnaser_map[np.argsort(-np.mean(hela_rnaser_map.T,axis=0))].T
k562_norm_map = k562_norm_map[np.argsort(-np.mean(k562_norm_map.T,axis=0))].T
k562_rnaser_map = k562_rnaser_map[np.argsort(-np.mean(k562_rnaser_map.T,axis=0))].T

hela_norm_map_filted = hela_norm_map[:,0:200].T
hela_rnaser_map_filted = hela_rnaser_map[:,0:200].T
k562_norm_map_filted = k562_norm_map[:,0:200].T
k562_rnaser_map_filted = k562_rnaser_map[:,0:200].T

#################
# Figure 4A-B   #
#################
# Figure 4 cross validation level
def draw_crossvalidation():
    # Set the axis and labels
    xLabel = softlist
    yLabel = xLabel.copy()
    yLabel.reverse()

    data = []
    n = len(xLabel)
    for i in range(n-1,-1,-1):
        temp = []
        for j in range(n):
            a, b = hela_rnaser[i][['chrom', 'start', 'end']], hela_rnaser[j][['chrom', 'start', 'end']]
            # k = len(list(set(a).intersection(set(b))))/len(a)
            k = len(pd.merge(a,b,on=['chrom', 'start', 'end'])) / len(a)
            temp.append(k)

        data.append(temp)

    data1 = []
    n = len(xLabel)
    for i in range(n - 1, -1, -1):
        temp = []
        for j in range(n):
            a, b = k562_rnaser[i][['chrom', 'start', 'end']], k562_rnaser[j][['chrom', 'start', 'end']]
            # k = len(list(set(a).intersection(set(b))))/len(a)
            k = len(pd.merge(a, b, on=['chrom', 'start', 'end'])) / len(a)
            temp.append(k)

        data1.append(temp)

    # Create a canvas
    fig = plt.figure(figsize=(15,12), dpi=300)
    # Define the canvas with 2*1 and plot it in the first position
    ax = fig.add_subplot(111)
    # Define the axis and labels
    ax.set_yticks(range(len(yLabel)))
    ax.set_yticklabels(yLabel,fontsize=25)
    ax.set_xticks(range(len(xLabel)))
    ax.set_xticklabels(xLabel,fontsize=25)
    ax.grid(False)
    # plot the clustermap by selecting specific colormap style
    im = ax.imshow(data, cmap=cmap1,vmin=0, vmax=1)
    # ax.set_aspect(5)
    # sns.clustermap(norm_map, method='ward', metric='euclidean', row_cluster=True, col_cluster=False, cmap=plt.cm.RdYlGn_r, \
    #         xticklabels=False, yticklabels=['find_circ','CIRCexplorer2','CircRNAfinder'], fontsize=15)
    for tick in ax.get_xticklabels():
        tick.set_rotation(60)
    # Add the right colorbar
    cb = plt.colorbar(im,shrink=0.5)
    cb.ax.tick_params(labelsize=20)
    font = {
        'family': ' sans-serif',
        'color': 'black',
        'weight': 'normal',
        'size': 20
    }
    cb.set_label('Proportion of common candidates', fontdict=font)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('fig3A.svg', dpi=600, bbox_inches='tight')
    plt.show()
    # Define the canvas with 2*1 and plot it in the first position
    fig = plt.figure(figsize=(15, 12), dpi=300)
    ax = fig.add_subplot(111)
    # Define the axis and labels
    ax.set_yticks(range(len(yLabel)))
    ax.set_yticklabels(yLabel, fontsize=25)
    ax.set_xticks(range(len(xLabel)))
    ax.set_xticklabels(xLabel, fontsize=25)
    ax.grid(False)
    # plot the clustermap by selecting specific colormap style
    im = ax.imshow(data1, cmap=cmap1, vmin=0, vmax=1)
    for tick in ax.get_xticklabels():
        tick.set_rotation(60)
    # Add the right colorbar
    cb1 = plt.colorbar(im,shrink=0.5)
    cb1.ax.tick_params(labelsize=20)

    cb1.set_label('Proportion of common candidates', fontdict=font)
    # show
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('fig3B.svg', dpi=600, bbox_inches='tight')
    plt.show()

figure_4 = draw_crossvalidation()


################
# Figure 5A-B   #
################
#figure 5A-B 
from matplotlib.colors import LinearSegmentedColormap
lut=dict(zip(np.unique(typelist),'gb'))
colcolors=[lut[i] for i in typelist]
cmap1=LinearSegmentedColormap.from_list('mycmap',['royalblue','white','darkred'])
def draw_reads_level():
    # Create a canvas
    fig = plt.figure(dpi=300)
    # Define the canvas with 1*1
    ax = fig.add_subplot(211)
    # Define the y-axis and labels
    ax.set_yticks(range(len(yLabel)))
    ax.set_yticklabels(yLabel)
    # ax.set_xticks(range(len(norm_map)))

    ax.grid(False)

    # plot the clustermap by selecting specific colormap style
    # im = ax.imshow(norm_map_filted, cmap=plt.cm.RdYlGn_r,aspect='auto')

    im = sns.clustermap(hela_norm_map_filted, row_cluster=False, col_colors=colcolors, xticklabels=softlist, yticklabels=False,\
                        cmap=cmap1, cbar_kws={'label': 'log2(read counts)','orientation': 'vertical','pad':0.15}, \
                        cbar_pos=(0.06, 0.2, 0.03, 0.4))
    plt.setp(im.ax_heatmap.xaxis.get_majorticklabels(),fontsize=20,rotation=45)
    sns.set('paper')
    # cb = plt.colorbar(im)
    # cb.ax.tick_params(labelsize=15)
    # font = {
    #     'family': ' sans-serif',
    #     'color': 'black',
    #     'weight': 'normal',
    #     'size': 15
    # }
    # cb.set_label('log2 (reads count)', fontdict=font)
    # plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('fig4A.svg', dpi=600,bbox_inches='tight')
    ax = fig.add_subplot(212)
    # Define the y-axis and labels
    ax.set_yticks(range(len(yLabel)))
    ax.set_yticklabels(yLabel)

    ax.grid(False)
    # im2 = ax.imshow(hypo_map_filted, cmap=plt.cm.RdYlGn_r,aspect='auto')
    im2 = sns.clustermap(hela_rnaser_map_filted, row_cluster=False, col_colors=colcolors, xticklabels=softlist, yticklabels=False,\
                        cmap=cmap1, cbar_kws={'label': 'log2(read counts)','orientation': 'vertical','pad':0.15}, \
                        cbar_pos=(0.06, 0.2, 0.03, 0.4))
    plt.setp(im2.ax_heatmap.xaxis.get_majorticklabels(),fontsize=20,rotation=45)
    # Add the right colorbar 
    # cb = plt.colorbar(im2)
    # cb.ax.tick_params(labelsize=15)
    # font = {
    #     'family': ' sans-serif',
    #     'color': 'black',
    #     'weight': 'normal',
    #     'size': 15
    # }
    # cb.set_label('log2 (reads count)', fontdict=font)
    # plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('fig4B.svg', dpi=600,bbox_inches='tight')
    plt.show()

figure_5_1 = draw_reads_level()

#################
# Figure 5C-D   #
#################
#figure 5C-D
def draw_reads_level_total():
    # Create a canvas
    fig = plt.figure(figsize=(20,20),dpi=300)
    # Define the canvas with 1*1
    ax = fig.add_subplot(211)
    # Set the y-axis and labels
    ax.set_yticks(range(len(yLabel)))
    ax.set_yticklabels(yLabel)
    # ax.set_xticks(range(len(norm_map)))

    ax.grid(False)
    # plot the clustermap by selecting specific colormap style
    # im = ax.imshow(norm_map, cmap=plt.cm.RdYlGn_r,aspect='auto')
    im = sns.clustermap(k562_norm_map_filted, row_cluster=False, col_colors=colcolors, xticklabels=softlist, yticklabels=False,\
                        cmap=cmap1, cbar_kws={'label': 'log2(read counts)','orientation': 'vertical','pad':0.15},\
                        cbar_pos=(0.06, 0.2, 0.03, 0.4))
    plt.setp(im.ax_heatmap.xaxis.get_majorticklabels(),fontsize=20,rotation=45)
    sns.set('paper')
    # cb = plt.colorbar(im)
    # cb.ax.tick_params(labelsize=15)
    # font = {
    #     'family': ' sans-serif',
    #     'color': 'black',
    #     'weight': 'normal',
    #     'size': 15
    # }
    # cb.set_label('log2 (reads count)', fontdict=font)
    ax = fig.add_subplot(212)
    # Set the y-axis and labels
    ax.set_yticks(range(len(yLabel)))
    ax.set_yticklabels(yLabel)
    # plt.gcf().subplots_adjust(bottom=0.15)

    plt.savefig('fig4C.svg', dpi=600, bbox_inches='tight')

    ax.grid(False)
    # im2 = ax.imshow(hypo_map, cmap=plt.cm.RdYlGn_r,aspect='auto')
    im2 = sns.clustermap(k562_rnaser_map_filted, row_cluster=False, col_colors=colcolors, xticklabels=softlist, yticklabels=False,\
                        cmap=cmap1, cbar_kws={'label': 'log2(read counts)','orientation': 'vertical','pad':0.15},\
                        cbar_pos=(0.06, 0.2, 0.03, 0.4))
    plt.setp(im2.ax_heatmap.xaxis.get_majorticklabels(),fontsize=20,rotation=45)
    # Add the right colorbar
    # cb = plt.colorbar(im2)
    # cb.ax.tick_params(labelsize=15)
    # font = {
    #     'family': ' sans-serif',
    #     'color': 'black',
    #     'weight': 'normal',
    #     'size': 15
    # }
    # cb.set_label('log2 (reads count)', fontdict=font)
    # plt.gcf().subplots_adjust(bottom=0.15)

    plt.savefig('fig4D.svg', dpi=600, bbox_inches='tight')
    plt.show()

figure_5_2 = draw_reads_level_total()



#####################################################################################
# Figure 6/S5 The frequencies and distribution of circRNAs after RNase R digestion  #
#####################################################################################
# RNaseR+/RNaseR-
xLabel=['hela','k562']
def draw_freq_rnaser(fig_name,input_norm,input_rnaser,xLabel,key):
    norm, rnaser = [],[]
    enrich_rate = pd.DataFrame()
    enrich_reads_rate = []
    for i in range(len(yLabel)):
        # hela cell
        if key == 'shared':
            total=[]
            if xLabel == 'hela':
                cross = pd.merge(hela_norm_inter,hela_rnaser_inter,on=['chrom','start','end'],how='outer')
                cross_inner = pd.merge(hela_norm_inter,hela_rnaser_inter,on=['chrom','start','end'],how='inner')
                cross = cross.drop_duplicates(subset=['chrom','start','end'], keep='first')
            else:
                cross = pd.merge(k562_norm_inter, k562_rnaser_inter, on=['chrom', 'start', 'end'], how='outer')
                cross_inner = pd.merge(k562_norm_inter, k562_rnaser_inter, on=['chrom', 'start', 'end'], how='inner')
                cross = cross.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
            for j in range(len(cross)):
                total.append('circRNA_'+str(j))
            cross['circRNA_ID']=total
            norm.append(pd.merge(input_norm[i],cross.loc[:, ['chrom','start','end','circRNA_ID']],on=['chrom','start','end'],how='left'))
            rnaser.append(pd.merge(input_rnaser[i],cross.loc[:, ['chrom','start','end','circRNA_ID']],on=['chrom','start','end'],how='left'))
            uneffected = set(list(pd.merge(cross_inner,cross,on=['chrom','start','end'],how='inner')['circRNA_ID']))
            uneffected_data = rnaser[i][rnaser[i]['circRNA_ID'].isin(list(uneffected))].loc[:,['circRNA_ID','n_reads']]
            uneffected_data['type']='uneffected'
            enrich = set(list(pd.merge(norm[i],cross.iloc[:,0:3],on=['chrom','start','end'],how='left')['circRNA_ID']))-uneffected
            enrich_data = rnaser[i][rnaser[i]['circRNA_ID'].isin(list(enrich))].loc[:,['circRNA_ID','n_reads']]
            enrich_data.columns = ['circRNA_ID', 'n_reads']
            enrich_data['type'] = 'enrich'
            depleted = set(list(pd.merge(rnaser[i],cross.iloc[:,0:3],on=['chrom','start','end'],how='left')['circRNA_ID']))-uneffected
            depleted_data = norm[i][norm[i]['circRNA_ID'].isin(list(depleted))].loc[:,['circRNA_ID','n_reads']]
            depleted_data['type'] = 'depleted'
        elif key == 'total':
            total = []
            cross = pd.merge(input_norm[i], input_rnaser[i], on=['chrom', 'start', 'end'], how='outer')
            cross_inner = pd.merge(input_norm[i], input_rnaser[i], on=['chrom', 'start', 'end'], how='inner')
            cross = cross.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
            for j in range(len(cross)):
                total.append('circRNA_' + str(j))
            cross['circRNA_ID'] = total
            norm.append(pd.merge(input_norm[i], cross.loc[:, ['chrom', 'start', 'end', 'circRNA_ID']],
                               on=['chrom', 'start', 'end'], how='left'))
            rnaser.append(pd.merge(input_rnaser[i], cross.loc[:, ['chrom', 'start', 'end', 'circRNA_ID']],
                                 on=['chrom', 'start', 'end'], how='left'))
            uneffected = set(
                list(pd.merge(cross_inner, cross, on=['chrom', 'start', 'end'], how='inner')['circRNA_ID']))
            uneffected_data = rnaser[i][rnaser[i]['circRNA_ID'].isin(list(uneffected))].loc[:,
                              ['circRNA_ID', 'n_reads']]
            uneffected_data['type'] = 'uneffected'
            enrich = set(list(pd.merge(norm[i], cross.iloc[:, 0:3], on=['chrom', 'start', 'end'], how='left')[
                                  'circRNA_ID'])) - uneffected
            enrich_data = rnaser[i][rnaser[i]['circRNA_ID'].isin(list(enrich))].loc[:, ['circRNA_ID', 'n_reads']]
            enrich_data.columns = ['circRNA_ID', 'n_reads']
            enrich_data['type'] = 'enrich'
            depleted = set(list(pd.merge(rnaser[i], cross.iloc[:, 0:3], on=['chrom', 'start', 'end'], how='left')[
                                    'circRNA_ID'])) - uneffected
            depleted_data = norm[i][norm[i]['circRNA_ID'].isin(list(depleted))].loc[:, ['circRNA_ID', 'n_reads']]
            depleted_data['type'] = 'depleted'

        # Combine the expression matrix from different datasets
        total_data = pd.concat([uneffected_data,enrich_data,depleted_data],axis=0)
        total_data = total_data.sort_values(by='n_reads',axis = 0,ascending = False).iloc[0:10000,:]
        total_data.columns = ['circRNA_ID','n_reads','type']
        total_data['n_reads'] = np.log2(total_data['n_reads'])
        enrich_reads_rate.append(total_data)

        enrich_rate.loc[yLabel[i],'uneffected']=[len(uneffected),len(enrich),len(depleted)][0]
        enrich_rate.loc[yLabel[i], 'enrich'] = [len(uneffected), len(enrich), len(depleted)][1]
        enrich_rate.loc[yLabel[i], 'depleted'] = [len(uneffected), len(enrich), len(depleted)][2]

    temp=enrich_rate.T
    # Calculate the fraction of different types (depletion, uneffected and enriched)
    percentages = np.zeros((10, 3))
    col_sum = np.sum(temp, axis=0)
    for i in range(temp.shape[0]):
        for j in range(len(temp.iloc[i,])):
            percentages[j, i] = temp.iloc[i, j] / col_sum[j] * 100

###################################################################
# Figure 6 The frequencies of circRNAs after RNase R digestion    #
###################################################################
    fig=plt.figure(figsize=(12,8),dpi=300)
    ax = fig.add_subplot(111)
    plt.style.use('ggplot')
    plt.bar(x = temp.columns.values,height=temp.loc['depleted',:],color='darkred',label='depleted', \
            tick_label=softlist)
    plt.bar(x=temp.columns.values,height=temp.loc['uneffected',:],color='gray',label='uneffected', \
            tick_label=softlist,bottom=list(temp.loc['depleted',:]))
    plt.bar(x=temp.columns.values,height=temp.loc['enrich',:],color='green',label='enrich', \
            tick_label=softlist,bottom=list(temp.loc['uneffected',:]+temp.loc['depleted',:]))
    plt.xticks(rotation=45)
    plt.legend(loc='best',fontsize=14)
    plt.tick_params(axis='both', labelsize=14)
    plt.ylabel('#circRNA',fontsize=20)
    plt.xlabel('software',fontsize=20)

    # Add the fraction into figure for stacked bar plot
    # search all of the bar segments and annotate
    bl = np.zeros(temp.shape[1])
    for j in [2,0,1]:

        for i,h in zip(range(len(temp.iloc[j,:])),temp.iloc[j,:]):
            x = i
            y = 0.5*h + bl[i]
            bl[i] = bl[i] + h
            ax.text(x, y, "{:.1f}%".format(percentages[i, j]), ha='center',fontsize=14)

    # plt.title(xLabel,fontsize=14)
    plt.savefig(fig_name+'_1.svg',dpi=600,bbox_inches='tight')
    plt.show()

################################################################################
# Figure S5 The enrichiment or depletion of circRNAs at reads level            #
################################################################################
    # Plot the distribution of circRNA ranked by expression level before and after RNaseR treat
    fig,ax=plt.subplots(figsize=(12,8),nrows=3,ncols=3,dpi=300)
    n=0
    for r in range(3):
        for c in range(3):
            ax[r,c].set(title=yLabel[n])
            data = enrich_reads_rate[n].reset_index(drop=True)
            ax[r, c].bar(data[data['type']=='enrich'].index,data[data['type']=='enrich']['n_reads'],linewidth=0,label='enrich',color='green')
            ax[r, c].bar(data[data['type'] == 'uneffected'].index, data[data['type'] == 'uneffected']['n_reads'],linewidth=0,
                         label='uneffected',color='gray')
            ax[r, c].bar(data[data['type'] == 'depleted'].index, data[data['type'] == 'depleted']['n_reads'],linewidth=0,
                         label='depleted',color='darkred')
            ax[r, c].set_ylabel('log2(read counts)',fontsize=9)
            x = (data.index.max())/2
            y = 10
            s = str(round(len(data[data['type'] == 'depleted'])*100/(len(data)),2))+'%'
            ax[r, c].text(x, y, s, fontdict=None,fontsize=15,ha='center',color='darkred')
            n += 1
    plt.savefig(fig_name+'_2.svg',dpi=300,bbox_inches='tight')
    plt.legend()
    plt.show()
    return enrich_reads_rate

figure_S5_1 = draw_freq_rnaser('fig_8A',hela_norm,hela_rnaser,xLabel[0],key='total')
figure_S5_2 = draw_freq_rnaser('fig_8B',k562_norm,k562_rnaser,xLabel[1],key='total')

figure_S5_3 = draw_freq_rnaser('fig_8C',hela_norm,hela_rnaser,xLabel[0],key='shared')
figure_S5_4 = draw_freq_rnaser('fig_8D',k562_norm,k562_rnaser,xLabel[1],key='shared')






