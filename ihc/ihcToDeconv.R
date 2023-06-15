##IHC to cell type

source("loadOrganoidData.R")

##first we download the data from the two tables
deconv<-sync$tableQuery("select * from syn23667404")$asDataFrame()
ihc<-sync$tableQuery("select * from syn24175711")$asDataFrame()


#join the data, should automatically work by specimenID
combined<-ihc%>%left_join(deconv)

##now compute correlation between both
cor.res<-combined%>%
  group_by(proteinMarker,cell_type,method)%>% ## for each marker
  summarize(numSamps=n(),corVal=cor(Percent_DAB_Stained,score,method='spearman'))%>%
  arrange(desc(corVal))

sigs<-subset(cor.res,abs(corVal)>0.65)


##some things are good! let's plot
library(ggplot2)
res=combined%>%
  subset(!is.na(method))%>%
  subset(cell_type%in%sigs$cell_type)%>%
  mutate(`Cell Type`=as.factor(cell_type))%>%
  group_by(method)%>%
  do(plots=
  ggplot(data=.,aes(y=score,x=Percent_DAB_Stained,col=`Cell Type`,shape=proteinMarker))+
  geom_point(aes(size=2))+
    scale_color_viridis_d()+
  facet_grid(proteinMarker~`Cell Type`)
    + ggtitle(unique(.$method)))


##then save files
apply(res,1,function(x) ggsave(paste0(x$method,'stainingCor.png'),x$plots))
