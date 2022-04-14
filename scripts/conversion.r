#!/usr/lib/R/bin/Rscript
############################################################################################
#  R script for convert the raw results into bed format on different datasets              #
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
	make_option(c("-r", "--remap"), type = "integer", default = FALSE,
	      action = "store", help = "The output file is the remapped and combined bed file! (only used for URBORUS)"),
	
	make_option(c("-a", "--algorithm"), type = "character", default = FALSE,
	      action = "store", help = "The input bed file is based on which algorithm.
		  FC: find_circ
		  CIRC2: CIRCexplorer2
		  CIRI2: CIRI2
		  CDBG: CircDBG
		  CF: CircRNAfinder
		  CM: CircMarker
		  MP: Mapsplice
		  KF: KNIFE
		  UB: URBORUS
		  DCC:DCC
		  SE:Segemehl"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# print(opt)

working_dir = opt$w
# get the directory contains raw prediction results of specific circRNA prediction software
data_dir = file.path(working_dir,opt$a)
# raw_data_file = file.path(data_dir,)

filenames=read.table(file.path(opt$i),sep='')

for (i in 1:length(filenames[,1])) {
	a=as.character(filenames[,1][i])
	print(a)
	if (opt$a=='FC') {
	ipath=paste(data_dir,a,'/circ_candidates.bed',sep='')} else if (opt$a=='CIRC2') {
	ipath=paste(data_dir,a,'/circ_candidates.bed',sep='')} else if (opt$a=='CF'){
	ipath=paste(data_dir,a,'/',a,'s_filteredJunctions.bed',sep='')} else if (opt$a=='MP'){
	ipath=paste(data_dir,a,'/fusions_raw.txt',sep='')} else if (opt$a=='CDBG'){
	ipath=paste(data_dir,a,'/Detection_Result/Brief_sum.txt',sep='')} else if (opt$a=='CM'){
	ipath=paste(data_dir,a,'/Detection_Result/Brief_sum.txt',sep='')} else if (opt$a=='CIRI2'){
	ipath=paste(data_dir,a,'/out.ciri',sep='')} else if (opt$a=='NCL'){
	ipath=paste(data_dir,a,'/NCLscan_out.result.info',sep='')} else if (opt$a=='UB'){
	ipath=paste(data_dir,a,'/circRNA_list.txt',sep='')} else if (opt$a=='DCC'){
	ipath=paste(data_dir,a,'/Alignment/',a,'Chimeric.out.junction.circRNA',sep='')} else if (opt$a=='SE'){
	ipath=paste(data_dir,a,'/',a,'.sum.bed',sep='')} else if (opt$a=='KF'){
	ipath=paste(data_dir,a,'/circRNA/glmReports/','background','_1__circJuncProbs.txt',sep='')}
	#ipath=paste(data_dir,a,'/circRNAreads/combinedReports/','background','_1__circJuncProbs.txt',sep='')}
	print(ipath)
	if (opt$r==1) {
		ipath=paste(data_dir,a,'/remapped_circ_candidates_convert_remapped.bed',sep='')
		data=read.table(ipath)
		data=data[!duplicated(data[,4]),]
		data1=read.table(paste(data_dir,a,'/circ_candidates_convert.bed',sep=''))
		data=data.frame(data[,1:3],data1[,4:8])
		write.table(data,paste(data_dir, a, '/circ_candidates_convert-1.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
	}
	if (opt$a=='CIRI2') {
		data=read.table(ipath, comment.char='') 
		data=data[-1,] } else if (opt$a=='CM') { 
		data=read.table(ipath,sep='')} else if (opt$a=='CDBG') { 
		data=read.table(ipath,sep='')} else if (opt$a=='KF') { 
		data=read.table(ipath,sep='\t', header=F);data = data[-1,];data$V5=data$V2} else {
		data=read.table(ipath,sep='\t')}
	if (length(data[,1])<1){
		
		write.table(data,paste(data_dir, a, '/circ_candidates_convert.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
		if (opt$a=='UB') {
		write.table(data[,1:4],paste(data_dir, a, '/circ_candidates_convert_remapped.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
		}
		next
	}
	V4=list()
	V1=list()
	V2=list()
	V3=list()
	if (opt$a=='CDBG') {
	  data$'reads'=data$V4 } 
	else if (opt$a=='CM') {
	  data$'reads'=data$V4 }
	for (j in 1:length(data$V4)) {
	  
		if (opt$a=='CDBG') {
		  V4[j]=paste('hsa_circRNA_', j, '/', as.character(data$'reads'[j]),sep="") } else if (opt$a=='CM') {
		  V4[j]=paste('hsa_circRNA_', j, '/', as.character(data$'reads'[j]),sep="") } else if (opt$a=='CIRI2') {
		  V4[j]=paste(as.character(data$V1[j]), '/', as.character(data$V5[j]),sep="") } else if (opt$a=='NCL') {
		  V4[j]=gsub('[.]','/',as.character(data$V1[j])) } else if (opt$a=='UB') {
		  V4[j]=paste('chr',as.character(data$V1[j]),':',as.character(data$V2[j]),'-',as.character(data$V3[j]),'/',as.character(data$V7[j]),sep="") } else if (opt$a=='DCC') {
		  V4[j]=paste(as.character(data$V1[j]),':',as.character(data$V2[j]),'-',as.character(data$V3[j]),'/',as.character(data$V5[j]),sep="") } else if (opt$a=='KF') {
		  V4[j]=paste('circRNA','/',data$V2[j],sep='') } else {
		  V4[j] = paste(as.character(data$V4[j]),'/',as.character(data$V5[j]),sep="") }
		if (opt$a=='MP') {
		  V1[j] = strsplit(as.character(data$V1[j]),split = '~',fixed = T)[[1]][1]
		  } else if (opt$a=='CIRI2') {
		  V1[j] = paste('chr',as.character(data$V2[j]),sep="") } else if (opt$a=='NCL') {
		  V1[j] = paste('chr',as.character(data$V2[j]),sep="") } else if (opt$a=='KF') {
                                  V1[j] = paste(as.character(strsplit(data$V1[j],'|',fixed=T)[[1]][1]),sep='')
		start_end = as.numeric(c(strsplit(strsplit(data$V1[j],'|',fixed=T)[[1]][3],':')[[1]][2],strsplit(strsplit(data$V1[j],'|',fixed=T)[[1]][2],':')[[1]][2]))
		V2[j] =  min(start_end)
		V3[j] =  max(start_end) } else {
		  V1[j] = paste('chr',as.character(data$V1[j]),sep="")
		  }
	}

	
	if (opt$a=='CIRI2') {
	  data$V2=unlist(V1) } else if (opt$a=='NCL') {
	    data$V2=unlist(V1) } else if (opt$a=='KF') {
	      
                      data$V2=unlist(V2);data$V3=unlist(V3)       
	      data$V1=unlist(V1) } else {
	      data$V1=unlist(V1) }
	

		if (opt$a=='CIRI2') {
			data=data[,-1]
			
			data$V3=as.numeric(data$V3)-1
			V6=data$V5
			data$V6=V6
			data$V5=unlist(V4) } else if (opt$a=='NCL') {
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
	  data=data[,c(2,3)]
	  data=t(apply(data,1,sort,decreasing=FALSE))
	  data=as.data.frame(data)
	  colnames(data)=c('V2','V3')
	  data$V1=unlist(V1)
	  data$V4=unlist(V4)
	  data=data.frame(data$V1,data$V2,data$V3,data$V4)
	  
	  data$data.V2=as.numeric(data$data.V2)-1
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
	}  else if (opt$a=='DCC') {
	  data$V2=data$V2-1
	}
	
	write.table(data,paste(data_dir, a, '/circ_candidates_convert.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
	if (opt$a=='UB') {
	write.table(data[,1:4],paste(data_dir, a, '/circ_candidates_convert_remapped.bed',sep=''), row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
	}
}




