##cross-data score
#' the goal of this script is to evaluate the correlations across all data modalities
#' to determine the similarity of the various growth metrics
#' 

##this file compares the correlation of samples to the primary across IHC, methylation, and rnaseq
#source("loadOrganoidData.R")
source('orgPlottingFunctions.R')
library(dplyr)
library(ggplot2)
library(readr)
##get metadata file
reticulate::use_virtualenv('r-reticulate')
syn<-reticulate::import('synapseclient')
sync<-syn$login()

#get methyl cor
tab <- read_csv(sync$get('syn51713151')$path)

tab$specimenID<-stringr::str_replace(tab$SampleID,'--','-')|>
  stringr::str_replace_all('-plus-',"+")|>
  stringr::str_replace('_.+','')|>
  stringr::str_replace_all('-M-*Mammo',' M')|>
  stringr::str_replace_all('-S--StemPro',' S')|>
  stringr::str_replace_all('-D--DMEM',' D')

annotes<-read_csv('rnaSeqProcessing/rnaseqorgids.csv')|>
  dplyr::select(-c(id,name))|>
  dplyr::distinct()%>%
  dplyr::mutate(specimenID=str_replace_all(specimenID,'NF0007-4-D$','NF0007-4 D'))%>%
  dplyr::mutate(specimenID=str_replace_all(specimenID,'NF0007-4-S$','NF0007-4 S'))%>%
  dplyr::mutate(specimenID=str_replace_all(specimenID,'NF0007-4-M$','NF0007-4 M'))%>%
  dplyr::mutate(experimentalCondition=str_replace_all(experimentalCondition,'Cytkines,Mammo','Cytokines,Mammo'))|>
  dplyr::ungroup()

annotes$Media<-sapply(annotes$experimentalCondition,function(x){
  ifelse(length(grep('Mammo',x))>0,'Mammo',
         ifelse(length(grep('StemPro',x))>0,'StemPro',
                ifelse(length(grep('DMEM',x))>0,'DMEM','Tumor')))})

annotes$extras<-sapply(annotes$experimentalCondition,function(x)
  ifelse(length(grep('Cytokines',x))>0,ifelse(length(grep("Forskoline",x))>0,'Cytokines,Forskolin','Cytokines'),'None'))

#join with metadata from file view
#file.metadata<-sync$tableQuery("SELECT name,specimenID,individualID FROM syn11601495 WHERE assay = 'bisulfiteSeq'")$asDataFrame()

#now we have proper metadata
#res<-file.metadata%>%
#  rowwise()%>%
#  mutate(SampleID=gsub('_R.*fastq.gz','',name))%>%
#  left_join(tab)


vars <- c('individualID','specimenID','Similarity','dataType','extras',
          'Media')

methyl.cor<-tab%>%
  left_join(annotes)|>
  mutate(Similarity=sqrt(linearR2))%>%
  distinct()%>%
  mutate(dataType='RRBS')%>%
  select(-linearR2)%>%
  select(vars)

#get rnaseq cor
rtab <- read_csv(sync$get('syn51713188')$path)|>
  dplyr::select(-Patient)|>
  rename(specimenID='altID')|>#,individualID='Patient')%>%
  left_join(annotes)|>
  mutate(dataType='rnaSeq')%>%
  mutate(specimenID=as.character(specimenID))%>%
  #left_join(tibble::rownames_to_column(annotes,'specimenID'),by='specimenID')%>%
  select(vars)


#get ihc cor
# #TODO
# itab <- sync$tableQuery("SELECT * from syn24988958")$asDataFrame()%>%
#   rename(specimenID='altID',individualID='Patient')%>%
#   dplyr::select(specimenID,Similarity)%>%
#   mutate(specimenID=as.character(specimenID))%>%
#   mutate(dataType='IHC')%>%
#   left_join(tibble::rownames_to_column(annotes,'specimenID'),by='specimenID')%>%
#   select(vars)
# 

ctab <- read_csv(sync$get('syn51713150')$path)|>
  dplyr::rename(Similarity='corVal')|>
  left_join(annotes)|>
  mutate(dataType='flow cytometry')%>%
  mutate(specimenID=as.character(specimenID))%>%
  #left_join(tibble::rownames_to_column(annotes,'specimenID'),by='specimenID')%>%
  select(vars)


full.tab<-rbind(methyl.cor,rtab,ctab)%>%
  drop_na()%>%
  mutate(Media=as.factor(Media))

full.tab%>%group_by(Media,dataType)%>%summarize(medSim=median(Similarity))

##statis for individual datasets
full.tab%>%lm(Similarity ~ Media, data=.)%>%summary()

##now plot
p<-full.tab%>%
  subset(Media!='Tumor')|>
  subset(individualID!='NF0002')|>
  ggplot(aes(x=extras,y=Similarity))+
  geom_boxplot(aes(alpha=0.8),outlier.shape=NA)+
  geom_jitter(aes(color=dataType,shape=dataType),size=3)+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(~Media)+
  theme_classic()
  


pp<-full.tab%>%
  ggplot(aes(x=Media,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=Media),outlier.shape=NA)+
  geom_jitter(aes(color=Media,shape=dataType))+
  scale_fill_manual(values=media_pal)+scale_color_manual(values=media_pal)+facet_grid(extras~.)+coord_flip()

ggsave('dataTypeCorrelation.pdf',p,width=12)
sync$store(syn$File('combinedCorrelation.pdf',parentId='syn11376065'))
# 
# ##now plot
# p2<-full.tab%>%
#   ggplot(aes(x=dataType,y=Similarity))+
#   geom_boxplot(aes(alpha=0.8,fill=dataType),outlier.shape=NA)+
#   #geom_jitter(aes(color=dataType,shape=extras))+
#   scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~.)
# 
# 
# 
# 
# ##now plot
# p3<-full.tab%>%
#   ggplot(aes(x=dataType,y=Similarity))+
#   geom_boxplot(aes(fill=dataType),outlier.shape=NA)+
#  #geom_jitter(aes(shape=individualID,col=dataType,alpha=0.1))+
#   scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~extras)+coord_flip()
# 
# ggsave('altCorrelation2.pdf',p3,height=6)
# # sync$store(syn$File('altCorrelation2.pdf',parentId='syn11376065'))
# 
# ##now plot
# p3<-full.tab%>%
#   subset(individualID!='NF0002')%>%
#   ggplot(aes(x=dataType,y=Similarity))+
#   geom_boxplot(aes(fill=dataType),outlier.shape=NA)+
#   #geom_jitter(aes(shape=individualID,col=dataType,alpha=0.1))+
#   scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~extras)+coord_flip()
# 
# ggsave('altCorrelation2_no0002.pdf',p3)
# sync$store(syn$File('altCorrelation2_no0002.pdf',parentId='syn11376065'))

##table with mean. 
full.tab%>%group_by(Media,extras)%>%summarize(meanSim=mean(Similarity),medSim=median(Similarity))%>%arrange(desc(meanSim))