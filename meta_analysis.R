##cross-data score
#' the goal of this script is to evaluate the correlations across all data modalities
#' to determine the similarity of the various growth metrics
#' 

source("loadOrganoidData.R")
library(dplyr)
library(ggplot2)

#get methyl cor
tab <- sync$tableQuery("SELECT * FROM syn24216591")$asDataFrame()

#join with metadata from file view
file.metadata<-sync$tableQuery("SELECT name,specimenID,individualID FROM syn11601495 WHERE assay = 'bisulfiteSeq'")$asDataFrame()

annotes<-sync$tableQuery("select * from syn24216672")$asDataFrame()


#now we have proper metadata
res<-file.metadata%>%
  rowwise()%>%
  mutate(SampleID=gsub('_R.*fastq.gz','',name))%>%
  left_join(tab)


nann<-annotes%>%
  mutate(extras=stringr::str_replace_all(experimentalCondition,"Mammo,*",""))%>%
  mutate(extras=stringr::str_replace_all(extras,"DMEM,*",""))%>%
  mutate(extras=stringr::str_replace_all(extras,"StemPro,*",""))%>%
  mutate(extras=stringr::str_replace_all(extras,',$',''))%>%
  mutate(extras=stringr::str_replace_all(extras,"^$","None"))

methyl.cor<-res%>%
  select(individualID,specimenID,linearR2)%>%
  mutate(Similarity=sqrt(linearR2))%>%
  distinct()%>%
  mutate(dataType='RRBS')%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'NF0007-4-M$','NF0007-4 M'))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'NF0008-1-M$','NF0008-1 M'))%>%
  select(-linearR2)%>%left_join(nann)

#get rnaseq cor
rtab <- sync$tableQuery("SELECT * FROM syn24828132")$asDataFrame()%>%
  rename(specimenID='altID',individualID='Patient')%>%
  dplyr::select(specimenID,Similarity)%>%
  mutate(dataType='rnaSeq')%>%
  left_join(nann,by='specimenID')



#get ihc cor
#TODO



##now merge into a single table
vars <- c('individualID','specimenID','Similarity','dataType',
          'Media','extras')

full.tab<-rbind(methyl.cor,rtab)%>%
  select(vars)%>%
  drop_na()%>%
  mutate(Media=as.factor(Media))

##now plot
p<-full.tab%>%
  ggplot(aes(x=Media,y=Similarity,fill=extras,shape=dataType,color=extras))+
  geom_boxplot(aes(alpha=0.8),outlier.shape=NA)+
  geom_jitter()+scale_fill_viridis_d()+scale_color_viridis_d()
ggsave('combinedCorrelation.png',p)
sync$store(syn$File('combinedCorrelation.png',parentId='syn11376065'))
