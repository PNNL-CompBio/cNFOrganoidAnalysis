## pull methylation data and plot correlations


source("../loadOrganoidData.R")
#get data from anela's table
tab <- sync$tableQuery("SELECT * FROM syn24216591")$asDataFrame()

#join with metadata from file view
file.metadata<-sync$tableQuery("SELECT name,specimenID,individualID FROM syn11601495 WHERE assay = 'bisulfiteSeq'")$asDataFrame()


#now we have proper metadata
res<-file.metadata%>%
  rowwise()%>%
  mutate(SampleID=gsub('_R.*fastq.gz','',name))%>%
  right_join(tab)%>%
  rename(altID='specimenID')%>%
  left_join(annotes)##join with annotations pulled from synapse

ddf<-res%>%
  select(individualID,altID,numCpGs10x,linearR2,numDMRs05,numHyperDMRs,CHmeth,Media,Cytokines,Forskoline)%>%
  distinct()

vars=c('linearR2','numDMRs05','numHyperDMRs','CHmeth')

library(cowplot)
library(ggplot2)
for(var in vars){
  
plist<-ddf%>%
  rename(Metric=var)%>%
  group_by(individualID)%>%
  group_map(.keep=TRUE,~ggplot(.x,aes(x=Media,y=Metric,shape=Forskoline,col=Cytokines))+
                  geom_point(aes(size=10))+scale_colour_manual(values=pal)+ggtitle(.x$individualID))
cowplot::plot_grid(plotlist=plist)   
fname=paste0(var,'MethylationAnalysis.pdf')
ggsave(filename=fname,width=12)
sync$store(syn$File(fname,parentId='syn24216715'))
}
#
 
 # ggsave(filename=paste0(prefix,'corPlots.pdf'),res,width=10)
#  return(ddf)