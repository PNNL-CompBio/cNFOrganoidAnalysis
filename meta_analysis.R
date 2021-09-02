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

pal = wesanderson::wes_palette('Darjeeling1')

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

##now merge into a single table
vars <- c('individualID','specimenID','Similarity','dataType',
          'Media','extras')

methyl.cor<-res%>%
  select(individualID,specimenID,linearR2)%>%
  mutate(Similarity=sqrt(linearR2))%>%
  distinct()%>%
  mutate(dataType='RRBS')%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'NF0007-4-M$','NF0007-4 M'))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'NF0008-1-M$','NF0008-1 M'))%>%
  select(-linearR2)%>%left_join(nann)%>%
  select(vars)

#get rnaseq cor
rtab <- sync$tableQuery("SELECT * FROM syn24828132")$asDataFrame()%>%
  rename(specimenID='altID',individualID='Patient')%>%
  dplyr::select(specimenID,Similarity)%>%
  mutate(dataType='rnaSeq')%>%
  left_join(nann,by='specimenID')%>%
  select(vars)



#get ihc cor
#TODO
itab <- sync$tableQuery("SELECT * from syn24988958")$asDataFrame()%>%
  rename(specimenID='altID',individualID='Patient')%>%
  dplyr::select(specimenID,Similarity)%>%
  mutate(dataType='IHC')%>%
  left_join(nann,by='specimenID')%>%
  select(vars)


ctab <- sync$tableQuery('SELECT * from syn25954974')$asDataFrame()%>%
  dplyr::select(specID,Similarity='corVal',specimenID)%>%
  mutate(dataType='flow cytometry')%>%
  left_join(nann,by='specimenID')%>%
  select(vars)


full.tab<-rbind(methyl.cor,rtab,itab,ctab)%>%
  drop_na()%>%
  mutate(Media=as.factor(Media))

##now plot
p<-full.tab%>%
  ggplot(aes(x=extras,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=extras),outlier.shape=NA)+
  geom_jitter(aes(color=extras,shape=dataType))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~.)+coord_flip()

ggsave('combinedCorrelation.pdf',p,width=6)
sync$store(syn$File('combinedCorrelation.pdf',parentId='syn11376065'))

##now plot
p2<-full.tab%>%
  ggplot(aes(x=dataType,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=dataType),outlier.shape=NA)+
  geom_jitter(aes(color=dataType,shape=extras))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~.)+coord_flip()

ggsave('altCorrelation.pdf',p2,width=6)
sync$store(syn$File('altCorrelation.pdf',parentId='syn11376065'))


##now plot
p3<-full.tab%>%
  ggplot(aes(x=dataType,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=dataType),outlier.shape=NA)+
  geom_jitter(aes(color=dataType,shape=individualID))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~extras)+coord_flip()

ggsave('altCorrelation2.pdf',p3)
sync$store(syn$File('altCorrelation2.pdf',parentId='syn11376065'))

##table with mean. 
full.tab%>%group_by(Media,extras)%>%summarize(meanSim=mean(Similarity),medSim=median(Similarity))%>%arrange(desc(meanSim))