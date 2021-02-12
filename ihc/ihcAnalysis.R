##compare IHC across samples

source("loadOrganoidData.R")

ihc<-sync$tableQuery("select * from syn24175711")$asDataFrame()


##it's unclear how to show that IHC data is changing between organoids and primary samples

##format as matrix
imat<-ihc%>%
  dplyr::select(specimenID,proteinMarker,Fraction_DAB_stained)%>%
  distinct()%>%
    subset(proteinMarker!='Pan-Cytokeritin')%>%
  pivot_wider(names_from=proteinMarker,values_from=Fraction_DAB_stained)%>%
  tibble::column_to_rownames('specimenID')


##select out annotations
iannotes<-ihc%>%
  dplyr::select(individualID,specimenID)%>%
  mutate(altID=specimenID)%>%distinct()%>%
  left_join(annotes)%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')%>%
  select(-altID)%>%
  rename(altID='individualID')
#  separate(specimenID,into=c('altID','extra'),sep=' ',remove=F)%>%
#  mutate(altID=stringr::str_replace(altID,'-$',''))%>%
#  dplyr::select(-extra)%>%
#  distinct()

res<-plotCorrelationBetweenSamps(t(imat),iannotes,prefix='IHC')

##TODO: store table on Synapse
write.table(res,file='tmp.csv',sep=',',row.names=F)
sync$store(syn$build_table('IHC-based sample correlations','syn11374354','tmp.csv'))

sync$store(syn$File('IHCcorplots.pdf',parentId='syn24226005'))

##format as matrix
imat<-ihc%>%
  dplyr::select(specimenID,proteinMarker,Fraction_DAB_stained)%>%
  distinct()%>%
  pivot_wider(names_from=proteinMarker,values_from=Fraction_DAB_stained,values_fill=-1)%>%
  tibble::column_to_rownames('specimenID')
library(pheatmap)
pheatmap(imat,cellwidth = 10,cellheight=10,annotation_row =select(iannotes,altID),clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',clustering_method = 'ward.D2',
         filename='ihcHeatmap.pdf')
sync$store(syn$File('ihcHeatmap.pdf',parentId='syn24226005'))