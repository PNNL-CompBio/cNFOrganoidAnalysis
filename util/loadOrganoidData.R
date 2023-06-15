
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
source("orgPlottingFunctions.R")
if(!require(nationalparkcolors)){
  devtools::install_github("katiejolly/nationalparkcolors")
  library(nationalparkcolors)
}
library(pheatmap)
library(ggplot2)
if(!require(wesanderson)){
  install.packages("wesanderson")
  library(wesanderson)
}
##get metadata file
#condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"

reticulate::use_virtualenv('r-reticulate')

# reticulate::use_condaenv(condaenv)
synapseclient<-reticulate::import('synapseclient')
sync<-synapseclient$login()
##get annotations
#didn't work on mac:
#tabres = sync$tableQuery('SELECT zScore,totalCounts,specimenID,Symbol,individualID,experimentalCondition FROM syn22878645')$asDataFrame()
#tabres2 = sync$tableQuery('SELECT zScore,totalCounts,specimenID,Symbol,individualID,experimentalCondition FROM syn21222341')$asDataFrame()
##instead download syn22878645 to org_rnaseq.csv and syn21222341 to pat_rnaseq.csv
tabres = readr::read_csv('org_rnaseq.csv')|>
  dplyr::select(zScore,totalCounts,specimenID,Symbol,individualID,experimentalCondition )
tabres2 = readr::read_csv('pat_rnaseq.csv')|>
  dplyr::select(zScore,totalCounts,specimenID,Symbol,individualID,experimentalCondition )

rnaseq<-rbind(tabres, tabres2)%>%
  distinct()%>%
  mutate(specimenID=str_replace_all(specimenID,'NF0007-4-D$','NF0007-4 D'))%>%
  mutate(specimenID=str_replace_all(specimenID,'NF0007-4-S$','NF0007-4 S'))%>%
    mutate(specimenID=str_replace_all(specimenID,'NF0007-4-M$','NF0007-4 M'))%>%
  mutate(experimentalCondition=str_replace_all(experimentalCondition,'Cytkines,Mammo','Cytokines,Mammo'))

##remove bad species from annnotations and RNA seq
badSpecs=c("NF0009-1 M","NF0009-1-M+C","NF0002-8-19-M+C+F","NF0002-8-19-D+C",
           "NF0002-8-19-D+C+F","NF0002-8-19-S+C")

##05-01-23 added another sample
#badSpecs <- c(badSpecs,'NF0009-1')

##get annotations first - also fails on M1
#annotes<-sync$tableQuery("select * from syn24216672")$asDataFrame()
annotes <- readr::read_csv("annotes.csv")
##keep track of organoid
orgs<-subset(annotes,experimentalCondition!='NaN')%>%select(specimenID)

#now create various annotations
biga<-annotes%>% 
  subset(specimenID%in%orgs$specimenID)%>%
  dplyr::select(-c(individualID,experimentalCondition))%>%
  mutate(val=1) %>% 
  pivot_wider(names_from = Media,values_from = val) %>% 
  mutate(across(-c(Cytokines,Forskoline,specimenID), ~replace_na(.x, 0))) %>%
  mutate(across(-c(Cytokines,Forskoline,specimenID), ~ifelse(.x==1, TRUE,FALSE)))

pannotes<-rnaseq%>%select(specimenID,individualID,experimentalCondition)%>%
  distinct()%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')

nannotes<-apply(annotes,2,as.character)%>%
  as.data.frame()%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')


### get drug data befmore we update the annotations
#drugData = sync$tableQuery('SELECT * FROM syn26145552')$asDataFrame()%>%
#  left_join(annotes,by='specimenID')


##then relabel things for future plotting
annotes<-annotes%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')%>%
  mutate(extras=stringr::str_replace_all(experimentalCondition,"Mammo,*",""))%>%
  mutate(extras=stringr::str_replace_all(extras,"DMEM,*",""))%>%
  mutate(extras=stringr::str_replace_all(extras,"StemPro,*",""))%>%
  mutate(extras=stringr::str_replace_all(extras,',$',''))%>%
  mutate(extras=stringr::str_replace_all(extras,"^$","None"))%>%
  mutate(Media=stringr::str_replace_all(Media,'None','Tumor'))

#names(patient_pal)<-unique(annotes$individualID)

##get RNAseq data
rnaseq<-rnaseq%>%subset(!specimenID%in%badSpecs)

#
pats=rownames(pannotes)[grep('patient',rownames(pannotes))]


library(edgeR)
countmat<-tidyr::pivot_wider(dplyr::select(rnaseq,totalCounts,specimenID,Symbol),names_from='specimenID',values_from='totalCounts',values_fill=list(totalCounts=0.0))|>tibble::column_to_rownames('Symbol')
dge <- DGEList(countmat)                        # DGEList object created from the count data
dge1<-cpm(dge) ##normalize counts to 1 million
dge2 <- calcNormFactors(dge1, method = "TMM") ##calculate size factors


#dgList <- estimateCommonDisp(countmat)
#dgList <- estimateTagwiseDisp(dgList)
norm_counts.table <- t(t(dge1)*(dge2))

mat<-rnaseq%>%
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  tibble::column_to_rownames('Symbol')%>%as.matrix()

vars<-apply(mat,1,var,na.rm=T)%>%
  sort(decreasing=T)
  
