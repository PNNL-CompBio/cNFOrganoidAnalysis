##download data for GEO upload

#this file queries the data and metadata and downloads the file

source('loadOrganoidData.R')

##get file list 
res = sync$tableQuery('SELECT * FROM syn11601495 WHERE ( ( "assay" = \'bisulfiteSeq\' OR "assay" = \'rnaSeq\' ) )')$asDataFrame()%>%
  dplyr::select(name,dataType,dataSubtype,readPair,id,specimenID,fileFormat)

res$fileType=stringr::str_split_fixed(res$name,pattern="\\.",n=2)[,2]


##first get processed files
##get md5sum for each
proc<-subset(res,dataSubtype=='processed')

##then get raw files
raw <- res%>%subset(dataSubtype!='processed')%>%
  tidyr::pivot_wider(id_cols=c('dataType','specimenID','fileType','fileFormat','dataSubtype'),
                     names_from='readPair')