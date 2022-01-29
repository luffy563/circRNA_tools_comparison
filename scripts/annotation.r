#!/usr/lib/R/bin/Rscript
# input files names
# Annotated by gene symbol
library('optparse')

option_list = list(
	make_option(c("-s", "--single"), type = "integer", default = FALSE,
	      action = "store", help = "The input bed file is output based on single-end seq or paired-end seq!"),
	make_option(c("-a", "--algorithm"), type = "character", default = FALSE,
	      action = "store", help = "The input bed file is based on which algorithm
		  FC: find_circ
		  CIRC2: CIRCexplorer2
		  CIRI2: CIRI2
		  CDBG: CircDBG
		  CM: CircMarker
		  MP: Mapsplice
		  KF: KNIFE
	    NCL:NCLscan
	    UB:UROBORUS"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# print(opt)

# which circRNAs prediction software
if (opt$a=='FC') {
	path='output/find_circ/'} else if (opt$a=='CIRC2') {
	path='output/CIRCexplorer2/'} else if (opt$a=='CDBG') {
	path='output/CircDBG/'} else if (opt$a=='CM') {
	path='output/CircMarker/'} else if (opt$a=='MP') {
	path='output/Mapsplice/'} else if (opt$a=='KF') {
	path='output/KNIFE/'} else if (opt$a=='CF') {
	path='output/CircRNAfinder/'} else if (opt$a=='CIRI2') {
	path='output/CIRI2/' } else if (opt$a=='NCL') {
	path='output/NCLscan/'} else if (opt$a=='UB') {
	  path='output/UROBORUS/'
	}
# single-end seq or paired-end seq
if (opt$s==1) {
	filenames=read.table('SRR_list.txt')}  else {
	filenames=read.table('SRR_list-1.txt')
}

ref_hsa=read.csv('ref_circ/remapped_hsa_hg19_hela_circRNA.csv',header = T)
ref_gene=read.table('ref_genome/gh38_ref.txt',sep='\t')
namses=c('chrom','start','end','name','score','strand','thickStart',
         'thichEnd','itemRgb','exonCount','exonSizes','exonOffsets','readNumber',
         'circType','geneName','isoformName','index','flankIntron')

ref_hsa_names=c('chrom','start','end','strand','circRNA_ID','genomic_length',
                'spliced_seq_lengthl','annotation',	'best_transcript','gene symbol')

ref_gene_names=c('n','isoformName','chrom','strand','start','end','thickstart',
                 'thickend','exonCount','exonStart','exonEnd','unknown','geneSymbol','1','2','3')
colnames(ref_hsa)=ref_hsa_names
colnames(ref_gene)=ref_gene_names
attach(ref_gene)
genename=data.frame(isoformName,geneSymbol)
attach(ref_hsa)

circID=data.frame(chrom,start,end,circRNA_ID)
genename=genename[!duplicated(genename$isoformName),]
for (i in 1:length(filenames[,1])) {
	a=as.character(filenames[,1][i])
	print(a)
	ipath=paste(path,a,'/circularRNA_known.txt',sep='')
	print(ipath)
	data=read.table(ipath,sep='\t')
	colnames(data)=namses
	data_merge=merge(data,genename,by='isoformName')
	data_merge_1=merge(data_merge,circID,by=c('chrom','start','end'))
	write.csv(data_merge_1,paste(path,a,'/circRNA_known_annotated.csv',sep=''),row.names=FALSE)
}



