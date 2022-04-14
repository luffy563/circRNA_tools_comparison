#!/bin/Rscript

library(UpSetR)
library(ggupset)
library(tidyverse)
library(reshape2)

#############
# Figure S2 #
#############

## Set the working directory
setwd("G:/LHF/pred_circ/output/")
# read the dataset for plotting upset plot
background_dataset<-read.csv('./BackgroundPlotdf.csv',header = T)
# head(background_dataset)
# get software list
soft<-colnames(background_dataset)[2:12]

## get the union set for each software
background_dataset_union<-background_dataset
background_dataset_union[,2:12][background_dataset_union[,2:12]>0]<-1
## get the intersection set for each software
background_dataset_intersect<-background_dataset
background_dataset_intersect[,2:12][background_dataset_intersect[,2:12]<3]<-0
background_dataset_intersect[,2:12][background_dataset_intersect[,2:12]>=3]<-1

##############
# Figure S2A #
##############
## upset plot of background_dataset for union set of software on datasets
upset(background_dataset_union,nsets = 11,point.size = 3,order.by = "freq",
      nintersects = 30,
      line.size = 1,text.scale=c(2, 1.5, 
                                 2, 1.5, 2, 1.5))

##############
# Figure S2B #
##############
## upset plot of background_dataset for intersection set of software on datasets
upset(background_dataset_intersect,nsets = 11,point.size = 3,
      line.size = 1,text.scale=c(2, 1.5, 
                                 2, 1.5, 2, 1.5))

## read the raw expression dataset 
BackgroundRaw<-read.table('./BackgroundRaw.csv',header = T,sep=',')
head(BackgroundRaw)
# replace the missing value by zero
BackgroundRaw$readcounts[is.na(BackgroundRaw$readcounts)]<-min(BackgroundRaw$readcounts,na.rm = T)
# replace the outliers by lower or upper value based on Capping method
q1 <- quantile(BackgroundRaw$readcounts, 0.01)        #get the value of 1% percentile  
q99 <- quantile(BackgroundRaw$readcounts, 0.99)       #get the value of 99% percentile  
BackgroundRaw[BackgroundRaw$readcounts < q1,]$readcounts <- q1  
BackgroundRaw[BackgroundRaw$readcounts > q99,]$readcounts <-q99  
# summary(BackgroundRaw$readcounts)

library(reshape2)
BackgroundLong=dcast(BackgroundRaw,datasets+software~circRNA_id,fun.aggregate = sum)
Backgroundwidth=melt(BackgroundLong,id.vars = c("datasets","software"),variable.name = "circRNA_id",value.name = 'readcounts')
Backgroundwidth$datasets<-unlist(lapply(Backgroundwidth$datasets,function(x){strsplit(as.character(x),"_",fixed = T)[[1]][2]}))
Backgroundwidth$datasets<-factor(Backgroundwidth$datasets,levels=seq(5,30,5))

# head(Backgroundwidth)

# filt out circRNA not identified by each software in any control datasets
union_circRNA_list<-list()
for (i in 1:length(soft)) {
  # i=1
  # i
  # i=1
  # soft
  # i=9
  circ_identified_df<-Backgroundwidth[Backgroundwidth$software==soft[i],]%>%group_by(software,circRNA_id)%>%
    summarise(n = sum(readcounts))
  # circ_identified<-circ_identified_df$circRNA_id[circ_identified_df$n>0]
  circ_identified<-background_dataset_intersect$circRNA_id[background_dataset_intersect[,soft[i]]==1]
  union_circRNA_list[[i]]<-circ_identified
  idx<-unlist(lapply(circ_identified,function(x){which(Backgroundwidth$circRNA_id==x)}))
  if (i==1){
    temp<-Backgroundwidth[idx,]
    BackgroundwidthFilt<-temp[temp$software==soft[i],]
  }else{
    temp<-Backgroundwidth[idx,]

    BackgroundwidthFilt<-rbind(BackgroundwidthFilt,temp[temp$software==soft[i],])
  }
  # Backgroundwidth[idx,]<-
}
library(ggplot2)
library(ggpubr)
library(ggsci)
##############
# Figure S2C #
##############
## The change of expression of circRNAs for each software on background datasets with increased simulated depths
ggplot(BackgroundwidthFilt,aes(x=datasets,y=readcounts,fill=datasets))+geom_boxplot()+
  # stat_compare_means(aes(x=datasets,y=readcounts))+
  scale_fill_npg()+
  guides(color=guide_legend())+
  facet_wrap(~software,scales = 'free_y')+
  xlab('simulated depths')+
  theme_bw()+theme(
    legend.position = 'none',
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25),
    axis.text = element_text(angle=45,hjust=1,size=18),
    strip.text = element_text(size=15)
  )


library(corrplot)
library(dplyr)
library(psych)
##############
# Figure S2D #
##############
## Correlation of software with different simulated depths
allCircRNAid<-unique(unlist(union_circRNA_list))
BackgroundLong$mean<-rowMeans(BackgroundLong[,allCircRNAid],na.rm = T)
BackgroundMatrix<-BackgroundLong[,c('datasets','software','mean')] %>%
  dcast(datasets~software)
BackgroundMatrix[is.na(BackgroundMatrix)]<-0
corMatrix<-cor(BackgroundMatrix[,2:length(BackgroundMatrix)])

# calculate the pearson correlation coefficient
res<-corr.test(corMatrix,adjust = 'BH')
# plot the pearson correlation heatmap
corrplot(res$r,p.mat = res$p,type = 'upper',col=COL2('RdYlBu', 200),order='hclust',
         sig.level = 0.05,insig = 'label_sig',diag = F,
         tl.col='black',tl.cex=1,pch.cex = 2)
