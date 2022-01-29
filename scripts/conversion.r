#!/usr/lib/Rscript
# input files names
library('optparse')

option_list = list(
	make_option(c("-s", "--single"), type = "integer", default = FALSE,
	      action = "store", help = "The input bed file is output based on single-end seq or paired-end seq!"),
	make_option(c("-r", "--remap"), type = "integer", default = FALSE,
	      action = "store", help = "The output file is the remapped and combined bed file! (only used for URBORUS)"),
	
	make_option(c("-a", "--algorithm"), type = "character", default = FALSE,
	      action = "store", help = "The input bed file is based on which algorithm.
		  FC: find_circ
		  CIRC2: CIRCexplorer2
		  CIRI2: CIRI2
		  CDBG: CircDBG
		  CM: CircMarker
		  MP: Mapsplice
		  KF: KNIFE
		  NCL: NCLscan
		  UB: URBORUS"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# print(opt)

# which circRNAs prediction software
if (opt$a=='FC') {
	path='output/find_circ/'} else if (opt$a=='CIRC2') {
	path='output/tophat_fusion/'
	} else if (opt$a=='CDBG') {
	path='output/CircDBG/'
} else if (opt$a=='CM') {
	path='output/CircMarker/'
} else if (opt$a=='MP') {
	path='output/Mapsplice/'
} else if (opt$a=='KF') {
	path='output/KNIFE/'
} else if (opt$a=='CF') {
  path='output/CircRNAfinder/'
} else if (opt$a=='MP') {
  path='output/Mapsplice/'
} else if (opt$a=='CIRI2') {
  path='output/CIRI2/'
} else if (opt$a=='NCL') {
  path='output/NCLscan/'
} else if (opt$a=='UB') {
  path='output/UROBORUS/'
}
# single-end seq or paired-end seq
if (opt$s==1) {
	filenames=read.table('SRR_list.txt')}  else {
	filenames=read.table('SRR_list-1.txt')
	}




for (i in 1:length(filenames[,1])) {
	a=as.character(filenames[,1][i])
	print(a)
	if (opt$a=='FC') {
	ipath=paste(path,a,'/circ_candidates.bed',sep='')} else if (opt$a=='CIRC2') {
	  ipath=paste(path,a,'/circ_candidates.bed',sep='')} else if (opt$a=='CF'){
	ipath=paste(path,a,a,'s_filteredJunctions.bed',sep='')} else if (opt$a=='MP'){
	ipath=paste(path,a,'/fusions_raw.txt',sep='')} else if (opt$a=='CDBG'){
	ipath=paste(path,a,'/Detection_Result/Brief_sum.txt',sep='')} else if (opt$a=='CM'){
	ipath=paste(path,a,'/Detection_Result/Brief_sum.txt',sep='')} else if (opt$a=='CIRI2'){
	ipath=paste(path,a,'/out.ciri',sep='')} else if (opt$a=='NCL'){
	ipath=paste(path,a,'/NCLscan_out.result.info',sep='')} else if (opt$a=='UB'){
	ipath=paste(path,a,'/circRNA_list.txt',sep='')}
	print(ipath)
	if (opt$r==1) {
		ipath=paste(path,a,'/remapped_circ_candidates_convert_remapped.bed',sep='')
		data=read.table(ipath)
		data1=read.table(paste(path,a,'/circ_candidates_convert.bed',sep=''))
		data=data.frame(data[,1:3],data1[,4:8])
		write.table(data,paste(path, a, '/circ_candidates_convert-1.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
	}
	if (opt$a=='CIRI2') {
		data=read.table(ipath, comment.char='') 
		data=data[-1,]} else {
		data=read.table(ipath)}
	V4=list()
	V1=list()
	if (opt$a=='CDBG') {
	  data$'reads'=data$V4 } else if (opt$a=='CM') {
	  data$'reads'=data$V4 }
	for (j in 1:length(data$V4)) {
	  
		if (opt$a=='CDBG') {
		  V4[j]=paste('hsa_circRNA_', j, '/', as.character(data$'reads'[j]),sep="") } else if (opt$a=='CM') {
		  V4[j]=paste('hsa_circRNA_', j, '/', as.character(data$'reads'[j]),sep="") } else if (opt$a=='CIRI2') {
		  V4[j]=paste(as.character(data$V1[j]), '/', as.character(data$V5[j]),sep="") } else if (opt$a=='NCL') {
		  V4[j]=gsub('[.]','/',as.character(data$V1[j])) } else if (opt$a=='UB') {
		  V4[j]=paste('chr',as.character(data$V1[j]),':',as.character(data$V2[j]),'-',as.character(data$V3[j]),'/',
		  as.character(data$V7[j]),sep="") } else {
		  V4[j] = paste(as.character(data$V4[j]),'/',as.character(data$V5[j]),sep="") }
		if (opt$a=='MP') {
		  V1[j] = strsplit(as.character(data$V1[j]),split = '~',fixed = T)[[1]][1]
		  } else if (opt$a=='CIRI2') {
		  V1[j] = paste('chr',as.character(data$V2[j]),sep="") } else if (opt$a=='NCL') {
		  V1[j] = paste('chr',as.character(data$V2[j]),sep="") } else {
		  V1[j] = paste('chr',as.character(data$V1[j]),sep="")
		  }
	}
	

	
	if (opt$a=='CIRI2') {
	  data$V2=unlist(V1) } else if (opt$a=='NCL') {
	    data$V2=unlist(V1) } else {
	      data$V1=unlist(V1) }
	
	if (opt$a!='CIRC2') {
		if (opt$a=='CIRI2') {
			data=data[,-1]
			
			data$V3=as.numeric(data$V3)-1
			V6=data$V5
			data$V6=V6
			data$V5=unlist(V4) }
		if (opt$a=='NCL') {
			data=data[,-1]
			
			V6=data$V6
			data=data[,-c(3,4,5,6,7)]
			
			data$V4=V6
			
			data=t(apply(data,1,sort,decreasing=FALSE))
			data=as.data.frame(data)
			colnames(data)=c('V2','V3','V1')
			data=data.frame(data$V1,data$V2,data$V3)

			data$data.V2=as.numeric(data$data.V2)-1
			data$V4=unlist(V4) } else {
			data$V4=unlist(V4) }
		if (opt$a=='MP') {
		  V3=data$V2
		  data$V2=data$V3
		  data$V3=V3
		} else if (opt$a=='CDBG') {
		  V3=data$V2
		  data$V2=data$V3
		  data$V3=V3
		  data$V2=data$V2-1
		} else if (opt$a=='CM') {
		  V3=data$V2
		  data$V2=data$V3
		  data$V3=V3
		  data$V2=data$V2-1
		}
	}


	write.table(data,paste(path, a, '/circ_candidates_convert.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
	if (opt$a=='UB') {
	write.table(data[,1:4],paste(path, a, '/circ_candidates_convert_remapped.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
	}
}




