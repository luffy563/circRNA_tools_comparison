#!/usr/lib/R/bin/Rscript
############################################################################################
#  R script for merged the annotated results by known hg38 circRNAs on different datasets  #
############################################################################################
## Date: 2022-4-13
## Author: Hongfei Liu
## Contact: lhf563@nwafu.edu.cn


library('optparse')

option_list = list(
	make_option(c("-w", "--working_dir"), type = "character", default = FALSE,
	      action = "store", help = "The working directory of this study"),
	make_option(c("-i", "--id"), type = "character", default = FALSE,
	      action = "store", help = "The list file contains accession id of datasets"),
	make_option(c("-a", "--algorithm"), type = "character", default = FALSE,
	      action = "store", help = "The input bed file is based on which algorithm
		FC: find_circ
		CIRC2: CIRCexplorer2
		CIRI2: CIRI2
		CDBG: CircDBG
		CF: CircRNAfinder
		CM: CircMarker
		MP: Mapsplice
		KF: KNIFE
		UB:UROBORUS
		DCC:DCC
		SE:Segemehl"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# print(opt)

working_dir = opt$w
# get the directory contains raw prediction results of specific circRNA prediction software
data_dir = file.path(working_dir,opt$a)

# which circRNAs prediction software

subpath='circularRNA_known.txt'
if (opt$a=='KF') {
	subpath='circ_candidates_convert_liftOver.bed'
} 

# get different accession id of datasets
filenames=read.table(file.path(opt$i),sep='')

# read the known circRNA information of hg38 version
ref_hsa=read.table('./ref_circ/combine_final_hg38.txt',header = T,sep='\t')
#ref_hsa=read.csv('./ref_circ/remapped_hsa_hg19_hela_circRNA.csv',header = T)

# read the known gene information of hg38 version
ref_gene=read.table('ref_genome/hg38_ref.txt',sep='\t')
if (subpath=='circularRNA_known.txt') {
namses=c('chrom','start','end','name','score','strand','thickStart',
         'thichEnd','itemRgb','exonCount','exonSizes','exonOffsets','readNumber',
         'circType','geneName','isoformName','index','flankIntron')
}else {
namses=c('chrom','start','end','name','readcounts')
}
#ref_hsa_names=c('circRNA_ID','chrom','start','end','strand','genomic_length',
#               'spliced_seq_lengthl','annotation',	'best_transcript','gene symbol')
ref_hsa_names=c('circRNA_ID','chrom','start','end','strand','genomic_length',
                'spliced_seq_lengthl','cell','annotation','feature','best_transcript','gene symbol','study')


ref_gene_names=c('n','isoformName','chrom','strand','start','end','thickstart',
                 'thickend','exonCount','exonStart','exonEnd','unknown','geneSymbol','1','2','3')
colnames(ref_hsa)=ref_hsa_names
colnames(ref_gene)=ref_gene_names
attach(ref_gene)
genename=data.frame(isoformName,geneSymbol)
attach(ref_hsa)

circID=data.frame(chrom,start,end,circRNA_ID)
genename=genename[!duplicated(genename$isoformName),]
if (opt$a=='KF') {
for (i in 1:length(filenames[,1])) {
	a=as.character(filenames[,1][i])
	print(a)
	ipath=paste(data_dir,a,'/',subpath,sep='')
	print(ipath)
	if (file.size(ipath)[1]==0){
	print(names)
	temp=data.frame(colnames=namses)
	data=t(temp)
	} else{
	data=read.table(ipath,sep='\t')
	print(data)
	}
	#tryCatch(data=read.table(ipath,sep='\t'),error=function(w){data=data.frame()})
	colnames(data)=namses
	#data_merge=merge(data,genename,by='isoformName')
	data_merge_1=merge(data,circID,by=c('chrom','start','end'))
	write.csv(data_merge_1,paste(data_dir,a,'/circRNA_known_annotated.csv',sep=''),row.names=FALSE,quote=FALSE)
}}else {
for (i in 1:length(filenames[,1])) {
	a=as.character(filenames[,1][i])
	print(a)
	ipath=paste(data_dir,a,'/',subpath,sep='')
	print(ipath)
	if (file.size(ipath)[1]==0){
	temp=data.frame(colnames=namses)
	data=t(temp)	
	} else{
	data=read.table(ipath,sep='\t')
	}
	#tryCatch(data=read.table(ipath,sep='\t'),error=function(w){data=data.frame()})
	colnames(data)=namses
	data_merge=merge(data,genename,by='isoformName')
	data_merge_1=merge(data_merge,circID,by=c('chrom','start','end'))
	write.csv(data_merge_1,paste(data_dir,a,'/circRNA_known_annotated.csv',sep=''),row.names=FALSE,quote=FALSE)
}
}



