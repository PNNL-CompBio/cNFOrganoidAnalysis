require(reticulate)
library(readxl)
library(stringr)

##get metadata file
condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"

  
#reticulate::use_condaenv(condaenv)
syn=reticulate::import('synapseclient')
sync=syn$login()
return(sync)
##get annotations
tabres = sync$tableQuery('SELECT id,name,assay,specimenID,individualID,experimentalCondition,readPair FROM syn11601495 WHERE ( ( "assay" = \'rnaSeq\' ) )')
annotes=tabres$asDataFrame()

rnaseq=annotes
#subset(annotes,assay=='rnaSeq')
#bfseq=subset(annotes,assay=='bisulfiteSeq')


##get file
fpath = sync$get('syn22783281')$path
patData = readxl::read_xlsx(fpath)[-c(1:2),]
names(patData)<-c('sample','specIDs')

newid<-data.frame(stringr::str_split_fixed(patData$specIDs,pattern='-',3))
names(newid)=c('indID','tumors','expCond')
newid$expCond=stringr::str_trim(newid$expCond)

newid<-data.frame(newid)%>%
  dplyr::mutate(newExpCond=dplyr::case_when(expCond==""~"None",
                                                       expCond%in%c('M(Mammo)','M (Mammo)')~'Mammo',
                                                       expCond=='D (DMEM)'~'DMEM',
                                                       expCond=='S (StemPro)'~'StemPro',
                                                       expCond=='S+C'~'Cytokines,StemPro',
                                                       expCond=='S+C+F'~'Cytokines,Forskoline,StemPro',
                                                       expCond=='M+C'~'Cytkines,Mammo',
                                                       expCond=='M+C+F'~'Cytokines,Forskoline,Mammo',
                                                       expCond=='D+C'~'Cytokines,DMEM',
                                                       expCond=='D+C+F'~'Cytokines,DMEM,Forskoline',
                                                       expCond%in%c("6","19")~"None"))

patData<-cbind(patData,newid)%>%dplyr::mutate(sample=stringr::str_trim(sample))


##first process the rnaseq
df=data.frame(name=rnaseq$name,
              do.call(rbind,lapply(stringr::str_split(rnaseq$name,pattern='_'),
                                   function(x) list(sample=x[[1]],rp=x[[4]]))))
df$sample<-unlist(df$sample)
ndf<-df%>%dplyr::left_join(rnaseq)%>%
  dplyr::left_join(patData,by='sample')%>%
  dplyr::select(id,name,specimenID='specIDs',
                individualID='indID',
                rp,
                experimentalCondition='newExpCond')%>%
  dplyr::mutate(readPair=stringr::str_replace(rp,"R",""))%>%dplyr::select(-rp)

rownames(ndf)<-rownames(rnaseq)
write.csv(ndf,file='annotes.csv')
write.csv(rnaseq,file='oldannotes.csv') ##merge manually and uplooad
#schem=sync$get('syn11601495')
#sync$store(syn$Table(schem,ndf,etag=schem$etag))
