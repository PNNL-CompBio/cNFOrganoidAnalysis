##compare IHC across samples

source("../loadOrganoidData.R")

ihc<-sync$tableQuery("select * from syn24175711")$asDataFrame()


##it's unclear how to show that IHC data is changing between organoids and primary samples

##format as matrix
imat<-ihc%>%
  dplyr::select(specimenID,proteinMarker,Fraction_DAB_stained)%>%
  distinct()%>%
    subset(!proteinMarker%in%c('Pan-Cytokeratin','Toluidine Blue'))%>%
  pivot_wider(names_from=proteinMarker,values_from=Fraction_DAB_stained,values_fn=list(Fraction_DAB_stained=mean))%>%
  tibble::column_to_rownames('specimenID')%>%as.matrix()


##select out annotations
iannotes<-ihc%>%
  dplyr::select(individualID,specimenID)%>% 
  mutate(altID=specimenID)%>%distinct()%>%
  left_join(tibble::rownames_to_column(annotes,'specimenID'),by=c('specimenID','individualID'))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')%>%
  select(-altID)
#  separate(specimenID,into=c('altID','extra'),sep=' ',remove=F)%>%
#  mutate(altID=stringr::str_replace(altID,'-$',''))%>%
#  dplyr::select(-extra)%>%
#  distinct()

res<-plotCorrelationBetweenSamps(t(imat),iannotes[rownames(imat),],prefix='IHC')


##now plot
p<-res%>%
  ggplot(aes(x=extras,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=extras),outlier.shape=NA)+
  geom_jitter(aes(color=extras,shape=Patient))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~.)+coord_flip()+ggtitle('IHC sample correlation')

ggsave('combinedIHCCorrelation.pdf',p,width=8)
sync$store(syn$File('combinedIHCCorrelation.pdf',parentId='syn11376065'))



##TODO: store table on Synapse
write.table(res,file='tmp.csv',sep=',',row.names=F)
cors<-sync$tableQuery('SELECT * FROM syn24988958')
sync$delete(cors)
sync$store(syn$build_table('IHC-based sample correlations','syn11374354','tmp.csv'))

sync$store(syn$File('IHCcorplots.pdf',parentId='syn24226005'))

##format as matrix
imat<-ihc%>%
  dplyr::select(specimenID,proteinMarker,Fraction_DAB_stained)%>%
  distinct()%>%
  pivot_wider(names_from=proteinMarker,values_from=Fraction_DAB_stained,values_fn=list(Fraction_DAB_stained=mean),values_fill=-1)%>%
  tibble::column_to_rownames('specimenID')
library(pheatmap)
pheatmap(imat,cellwidth = 10,cellheight=10,annotation_row =select(iannotes,individualID),clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',clustering_method = 'ward.D2',
         filename='ihcHeatmap.pdf')
sync$store(syn$File('ihcHeatmap.pdf',parentId='syn24226005'))