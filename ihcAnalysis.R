##compare IHC across samples

source("loadOrganoidData.R")

ihc<-sync$tableQuery("select * from syn24175711")$asDataFrame()


##it's unclear how to show that IHC data is changing between organoids and primary samples

##format as matrix
imat<-ihc%>%
  dplyr::select(specimenID,proteinMarker,Fraction_DAB_stained)%>%
  distinct()%>%
  pivot_wider(names_from=proteinMarker,values_from=Fraction_DAB_stained)%>%
  tibble::column_to_rownames('specimenID')


##select out annotations
iannotes<-ihc%>%
  dplyr::select(individualID,specimenID)%>%
  mutate(altID=specimenID)%>%distinct()%>%
  left_join(annotes)%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')
#  separate(specimenID,into=c('altID','extra'),sep=' ',remove=F)%>%
#  mutate(altID=stringr::str_replace(altID,'-$',''))%>%
#  dplyr::select(-extra)%>%
#  distinct()
res<-plotCorrelationBetweenSamps(imat,iannotes,prefix='IHC')
