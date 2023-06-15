#Sys.setenv(RETICULATE_PYTHON = '/Users/bade228/opt/anaconda3/envs/r2/bin/python3')
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(nationalparkcolors)
library(pheatmap)
library(ggplot2)
library(magrittr)
library(Biobase)
library(GSVA)
library(reticulate)
library(stringi)
library(cowplot)
library(wesanderson)
##get metadata file
#condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"

#pal = c(wesanderson::wes_palette('Darjeeling1'),wesanderson::wes_palette('Darjeeling2'))

#syn=reticulate::
syn=import('synapseclient')
sync=syn$login()
#pal = nationalparkcolors::park_palette('Acadia')

 query_tables <- function(table1synID,table2synID,query1cols,query2cols) {
   mat_dat <- read.csv(sync$tableQuery(paste0('SELECT ',query1cols,' FROM ',table1synID))$filepath,sep=',',header=TRUE) %>%
              select('Gene_Set_Type','Description','Gene')
   mat_dat = mat_dat%>%rename(Symbol=Gene)
   rnaseq <- read.csv(sync$tableQuery(paste0('SELECT ',query2cols,' FROM ',table2synID))$filepath,sep=',',header=TRUE) %>%
             select('Symbol','zScore','specimenID','individualID','experimentalCondition','diagnosis')
 
   #table3 = rnaseq[which(rnaseq$experimentalCondition=='None'),] #&rnaseq$diagnosis=='Neurofibromatosis 1'
   combined<-left_join(mat_dat,rnaseq,by=c("Symbol"))
   return(combined)
 }
 combo.mat.genecounts <- query_tables('syn24183888','syn22878645','Gene_Set_Type,Description,Gene','Symbol,zScore,specimenID,individualID,experimentalCondition,diagnosis')
 
 names <- unique(combo.mat.genecounts$Gene_Set_Type)
# restruc <- combo.mat.genecounts %>% select(specimenID,Symbol,zScore) %>%
#            distinct(specimenID,Symbol,.keep_all=TRUE) %>% 
#            drop_na(specimenID) %>%
#            spread(specimenID,zScore)
# rownames(restruc) <- restruc[,1]
# restruc2 <- restruc[,2:ncol(restruc)]
# mat.rest <- data.matrix(restruc2,rownames.force=TRUE)

# +
source('loadOrganoidData.R')


gset <- read.table(sync$get('syn26199051')$path,header=TRUE,sep='\t')

mat.rest <- mat

GSET <- list()
for (x in names) {
  val = gset[[x]]
  GSET[[x]] = stri_remove_empty(val)
}

# -

res <- gsva(mat.rest,GSET,method="ssgsea",min.sz=1,max.sz=Inf,mx.diff=TRUE,parallel.sz=1,verbose=TRUE)
res <- t(as.data.frame(res))
print(res)

master_table <- as.data.frame(res) %>%
                tibble::rownames_to_column(var='specimenID') %>%
             #   mutate(specimenID=str_replace_all(specimenID,fixed('-D'),' D')) %>%
            #    mutate(specimenID=str_replace_all(specimenID,fixed('-S'),' S')) %>%
            #    mutate(specimenID=str_replace_all(specimenID,fixed('-M'),' M')) %>%
                set_rownames(.$specimenID)

colnames(master_table)<-stringr::str_replace_all(colnames(master_table),'NABA_','')
modtab <- master_table %>%
          mutate(experimentalCondition = case_when(grepl("M\\+C\\+F",specimenID) ~ "Cytokines,Forskoline,Mammo",
                                                   grepl("S\\+C\\+F",specimenID) ~ "Cytokines,Forskoline,StemPro",
                                                   grepl("D\\+C\\+F",specimenID) ~ "Cytokines,Forskoline,DMEM",
                                                   grepl("M\\+C",specimenID) ~ "Cytokines,Mammo",
                                                   grepl("S\\+C",specimenID) ~ "Cytokines,StemPro",
                                                   grepl("D\\+C",specimenID) ~ "Cytokines,DMEM",
                                                   grepl("M",specimenID) ~ "Mammo",
                                                   grepl("S",specimenID) ~ "StemPro",
                                                   grepl("D",specimenID) ~ "DMEM",
                                                   TRUE ~ "None"
                                                   )
          )
modtab <- modtab %>% 
          mutate(individualID=sapply(str_split(specimenID, " "), function(x) x[1]))

## add in heatmap here
patindiv<-nannotes[pats,'individualID']

###plot heatmap of correlations
sannotes<-nannotes%>%as.data.frame(stringsAsFactors=FALSE)%>%
  dplyr::mutate(cohort=ifelse(individualID%in%patindiv,'cNF','Organoid'))%>%
  dplyr::select(Media,Cytokines,Forskoline,individualID)%>%
  mutate(Cytokines=as.character(Cytokines))%>%
  mutate(Forskoline =as.character(Forskoline))%>%
  mutate(Media=stringr::str_replace_all(Media,'None','Tumor'))%>%
  mutate(Media=as.factor(Media))

rownames(sannotes)<-rownames(nannotes)

annote.colors<-lapply(names(sannotes),function(x) c(`FALSE`='white',`TRUE`='grey'))
names(annote.colors)<-names(sannotes)
annote.colors$Media<-media_pal#[1:length(unique(sannotes$Media))]#c(None='white',StemPro='black',Mammo='darkgrey',DMEM='lightgrey')
#nnote.colors$cohort<-c(cNF='white',Organoid='black')
annote.colors$individualID=patient_pal#c(patient_pal,pal)[1:length(unique(sannotes$individualID))]
#                             NF0009='lightgrey',NF0007='darkslategrey')
#names(annote.colors$individualID)<-unique(sannotes$individualID)


modtab%>%
  subset(!specimenID%in%pats)%>%
  dplyr::select(-c(specimenID,experimentalCondition,individualID,MATRISOME))%>%
  t()%>%
  as.matrix()%>%
  pheatmap(.,annotation_col = sannotes,annotation_colors=annote.colors,
           cellheight=10,cellwidth=10,filename='matrisomeGSVA.pdf')



modtab%>%
  dplyr::select(-c(specimenID,experimentalCondition,individualID,MATRISOME))%>%
  t()%>%
  as.matrix()%>%
  pheatmap(.,annotation_col = sannotes,annotation_colors=annote.colors,
           cellheight=10,cellwidth=10,filename='all_matrisomeGSVA.pdf')

mat <- master_table[,-1]
genelist <- colnames(mat)

full.tab<-do.call(rbind,lapply(genelist,function(group) {
  print(group)
    plotCorrelationBetweenSamps(mat=t(mat),annotes=annotes,prefix=group)
}))


##Now add in the correlation boxplots
##now plot
p<-full.tab%>%
  ggplot(aes(x=extras,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=extras),outlier.shape=NA)+
  geom_jitter(aes(color=extras,shape=dataType))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~.)+coord_flip()

ggsave('combinedMatrisomeCorrelation.pdf',p,width=8)
sync$store(syn$File('combinedMatrisomeCorrelation.pdf',parentId='syn11376065'))

##now plot
p2<-full.tab%>%
  ggplot(aes(x=dataType,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=dataType),outlier.shape=NA)+
  geom_jitter(aes(color=dataType,shape=extras))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~.)+coord_flip()

ggsave('altMatrisomeCorrelation.pdf',p2,width=8)
sync$store(syn$File('altMatrisomeCorrelation.pdf',parentId='syn11376065'))


##now plot
p3<-full.tab%>%
  ggplot(aes(x=dataType,y=Similarity))+
  geom_boxplot(aes(alpha=0.8,fill=dataType),outlier.shape=NA)+
  geom_jitter(aes(color=dataType,shape=Patient))+
  scale_fill_manual(values=pal)+scale_color_manual(values=pal)+facet_grid(Media~extras)+coord_flip()

ggsave('altMatrisomeCorrelation2.pdf',p3,width=12)
sync$store(syn$File('altMatrisomeCorrelation2.pdf',parentId='syn11376065'))

