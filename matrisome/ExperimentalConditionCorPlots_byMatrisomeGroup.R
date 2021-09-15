Sys.setenv(RETICULATE_PYTHON = '/Users/bade228/opt/anaconda3/envs/r2/bin/python3')
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
syn=reticulate::import('synapseclient')
sync=syn$login()
pal = nationalparkcolors::park_palette('Acadia')

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
restruc <- combo.mat.genecounts %>% select(specimenID,Symbol,zScore) %>%
           distinct(specimenID,Symbol,.keep_all=TRUE) %>% 
           drop_na(specimenID) %>%
           spread(specimenID,zScore)
rownames(restruc) <- restruc[,1]
restruc2 <- restruc[,2:ncol(restruc)]
mat.rest <- data.matrix(restruc2,rownames.force=TRUE)

gset <- read.table(sync$get('syn26199051')$path,header=TRUE,sep='\t')
gset <- gset[2:nrow(gset),]
names <- colnames(gset)
GSET <- list()
for (x in names) {
  val = gset[[x]]
  GSET[[x]] = stri_remove_empty(val)
}

res <- gsva(mat.rest,GSET,method="gsva",min.sz=1,max.sz=Inf,mx.diff=TRUE,parallel.sz=1,verbose=TRUE)
res <- t(as.data.frame(res))
print(res)

master_table <- as.data.frame(res) %>%
                tibble::rownames_to_column(var='specimenID') %>%
                mutate(specimenID=str_replace_all(specimenID,fixed('-D'),' D')) %>%
                mutate(specimenID=str_replace_all(specimenID,fixed('-S'),' S')) %>%
                mutate(specimenID=str_replace_all(specimenID,fixed('-M'),' M')) %>%
                set_rownames(.$specimenID)

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

annotes<-modtab%>%
         dplyr::select(specimenID,individualID,experimentalCondition)%>%
         separate(specimenID,into=c('altID','extra'),sep=' ',remove=F,fill='right')%>%
         mutate(altID=stringr::str_replace(altID,'-$',''))%>%
         dplyr::select(-extra)%>%
         distinct() %>%
         drop_na(specimenID) %>%
         tibble::remove_rownames() %>%
         tibble::column_to_rownames(var='specimenID')
annotes$Media<-rep("None",nrow(annotes))
annotes$Media[grep("DMEM",annotes$experimentalCondition)]<-'DMEM'
annotes$Media[grep("StemPro",annotes$experimentalCondition)]<-'StemPro'
annotes$Media[grep('Mammo',annotes$experimentalCondition)]<-'Mammo'

pats=annotes$altID[grep('patient',annotes$individualID)]

specs<-c(rownames(annotes)[grep('patient',rownames(annotes))],
         'NF0002-8-19 M','NF0009-1- M+C+F','NF0012-3-6 M')

for(x in c('kines','Forskoline')){
  annotes[[x]]<-FALSE
  annotes[[x]][grep(x,annotes$experimentalCondition)]<-TRUE
  # annotes[[x]]=as.character(annotes)
}

annotes<-annotes%>%dplyr::rename(Cytokines='kines')#%>%dplyr::select(-experimentalCondition)

mat <- master_table[,1:10]
genelist <- colnames(mat)
my_geneset <- unique(modtab$individualID)
alldat <- cbind(annotes,mat)

plotCorrelationBetweenSamps<-function(mat,annotes,prefix){
  # filters out patient names
  orgs<-setdiff(annotes$altID,pats)
  # applied to every individualID in list
  dlist <- lapply(my_geneset,function(x) {
    alldat = alldat[which(alldat$individualID==x),]
    alldat = alldat %>% dplyr::select(Media,Cytokines,Forskoline,prefix)
    colnames(alldat)[4] <- 'Similarity'
    return(alldat)
  })
  
  #write.csv(dlist,paste0(prefix,'orgCorrelations.csv'),row.names=F)
  print(dlist)
  names(dlist)<-orgs
  plist<-lapply(orgs,function(pat){
    pdat<-dlist[[pat]]
    
    names(pdat)[names(pdat)==pat] = "Similarity"
    
    pdat%>%ggplot(aes(y=Similarity,x=Media,shape=Forskoline,color=Cytokines))+
      geom_point(aes(size=10))+scale_colour_manual(values=pal)+
      ggtitle(pat)
  })
  
  
  library(cowplot)
  res=cowplot::plot_grid(plotlist=plist)
  ggsave(filename=paste0('experimentalConditionCorPlots/',prefix,'corPlots.pdf'),res,width=10)
  #return(ddf)
}

for (group in genelist) {
    plotCorrelationBetweenSamps(mat=mat,annotes=annotes,prefix=group)
}



