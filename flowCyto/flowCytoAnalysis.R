##compare flowCytometry across samples

source("../loadOrganoidData.R")
#library(tidyverse)
#ihc<-sync$tableQuery("select * from syn24175711")$asDataFrame()
flow <- sync$get('syn25953258')$path%>%read.csv2(sep=',',skip=1,header=F)
colnames(flow)<-c('specimenID','proteinMarker','well1','well2','well3')

flow<-flow%>%
  #rowwise()%>%
  #mutate(meanVal=mean(as.numeric(c(val1,val2,val3)),na.rm=T))%>%
  subset(specimenID!="")%>%
  pivot_longer(c(well1,well2,well3),names_to='replicate',values_to='value')%>%
  mutate(specID=paste(specimenID,replicate,sep='_'))%>%
  mutate(value=as.numeric(value))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'-M$',' M'))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'-D$',' D'))%>%
           mutate(specimenID=stringr::str_replace_all(specimenID,'-S$',' S'))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'NF00012','NF0012'))
  
##it's unclear how to show that IHC data is changing between organoids and primary samples

##format as matrix
fmat<-flow%>%
  dplyr::select(specID,proteinMarker,value)%>%
  replace_na(list(value=0.0))%>%
  distinct()%>%
  pivot_wider(names_from=proteinMarker,values_from=value,values_fn=list(value=mean),
              values_fill=list(value=0.0))%>%
  tibble::column_to_rownames('specID')


##compute correlation
corvals <- cor(t(fmat))%>%as.data.frame()%>%
  tibble::rownames_to_column('origSpecId')%>%
  pivot_longer(cols=c(-origSpecId),names_to='specID',values_to='corVal')

fannotes<-sync$tableQuery('select * from syn24216672')$asDataFrame()%>%
  #tibble::remove_rownames()%>%
  #tibble::column_to_rownames('specimenID')#%>%
  #select(-experimentalCondition)
  right_join(flow)%>%
  select(individualID,experimentalCondition,Media,Cytokines,Forskoline,specimenID,specID,replicate)%>%
  distinct()

toplot<-corvals%>%##TODO: filter for equal replicates
  left_join(fannotes)%>%
  subset(experimentalCondition=='None')%>%
  rename(firstRep='replicate',firstIndiv='individualID')%>%
  select(-c(experimentalCondition,specID,Media,Cytokines,Forskoline,specimenID))%>%
  rename(specID='origSpecId')%>%left_join(fannotes)%>%
  subset(experimentalCondition!='None')%>%distinct()%>%
  mutate(equal=(firstRep==replicate&firstIndiv==individualID))%>%
  subset(equal)%>%
  subset(!is.na(corVal))

toplot%>%distinct()%>%
  ggplot(aes(x=Media,y=corVal,col=Cytokines,shape=Forskoline,size=5))+geom_point()+scale_colour_manual(values=pal)+facet_grid(~individualID)
#  distinct()

ggsave('flowCytocorplots.pdf')

##TODO: store table on Synapse
write.table(toplot,file='tmp.csv',sep=',',row.names=F)
cors<-sync$tableQuery('SELECT * FROM syn25954974')
sync$delete(cors)
sync$store(syn$build_table('Flow Cytometry-based sample correlations','syn11374354','tmp.csv'))

sync$store(syn$File('flowCytocorplots.pdf',parentId='syn24226005'))

mannotes<-fannotes%>%tibble::column_to_rownames('specID')
library(pheatmap)
pheatmap(fmat,cellwidth = 10,cellheight=10,
         annotation_row =select(mannotes,replicate,individualID),
         clustering_method = 'ward.D2',
         filename='flowCytoHeatmap.pdf')
#sync$store(syn$File('flowCytoHeatmap.pdf',parentId='syn24226005'))