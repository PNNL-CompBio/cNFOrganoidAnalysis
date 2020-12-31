##run tumor deconvolution on organoids


library(immunedeconv)  
source("loadOrganoidData.R")
library(ggplot2)

#' runImmuneDeconv
#' @param tab
#' @param method
#' @prefix 
#' @return Tidied data frame
runImmuneDeconv<-function(tab,method,prefix='organoids'){
  #run MCP counter
  mat<-reshape2::acast(tab,Symbol~specimenID,value.var='totalCounts',fun.aggregate=mean,na.rm=T)
  nas<-which(apply(mat,1,function(x) any(is.na(x))))
  if(length(nas)>0)
    mat<-mat[-nas,]
  res<-deconvolute(mat,method)
  
  df<-dplyr::select(tab,c(specimenID,experimentalCondition))%>%
      unique()#%>%
  #    rename(study='studyName')
  rownames(df)<-df$specimenID
  #save as heatmap with metadata
  mtab<-res%>%select(-cell_type)%>%as.data.frame()
  rownames(mtab)<-res$cell_type
  library(pheatmap)

  ##now tidy up data to table
  td<-tidyr::gather(res,key="specimenID",value="score",-cell_type )%>%
    left_join(df,by='specimenID')
  td$method=method
  return(td)
}
##NEED TO DOWNLOD FILES MANUALLY I CANNOT DISTRIBUTE
set_cibersort_binary('./CIBERSORT.R')
set_cibersort_mat('./LM22.txt')

xc<-runImmuneDeconv(rnaseq,'xcell')
mc<-runImmuneDeconv(rnaseq,'mcp_counter')
cs<-runImmuneDeconv(rnaseq,'cibersort')


full.tab<-rbind(xc,mc,cs)

write.table(full.tab,'file.tsv',sep='\t',quote=FALSE,row.names=F)
tab<-syn$build_table('cNF Organoid Deconvolution','syn11374354','file.tsv')
sync$store(tab)

condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
plot.annotes<-biga%>%
    select(c(Cytokines,Forskoline,DMEM,Mammo,StemPro,specimenID))%>%
    mutate(Cytokines=as.character(Cytokines))%>%
  mutate(Forskoline=as.character(Forskoline))%>%
  mutate(DMEM=as.character(DMEM))%>%mutate(Mammo=as.character(Mammo))%>%mutate(StemPro=as.character(StemPro))%>%
  tibble::column_to_rownames('specimenID')

annote.colors<-lapply(names(plot.annotes),function(x) c(`FALSE`='white',`TRUE`='black'))
names(annote.colors)<-names(plot.annotes)


##MCP Counter
mat<-mc%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
    tibble::column_to_rownames('cell_type')
  

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,
         cellwidth=10,annotation_colors=annote.colors,filename = 'mcpTypes.pdf',width=18)

all.d<-lapply(condlist,function(x){
  getDifferencesInCondition(mat,biga,x,'mcpTypes',doPlot=FALSE)  
})
names(all.d)<-condlist
cellCors<-data.frame(all.d)%>%mutate(dataType='mcpCounterPreds')
ddf<-plotCorrelationBetweenSamps(mat,annotes,'mcpCounter')

##CIBERSORT

mat<-cs%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

all.c<-lapply(condlist,function(x){
    getDifferencesInCondition(mat,biga,x,'Cibersort',doPlot=FALSE)  
})
names(all.c)<-condlist
cellCors<-data.frame(all.c)%>%mutate(dataType='cibersortPreds')%>%rbind(cellCors)
ddf<-plotCorrelationBetweenSamps(mat,annotes,'cibersort')

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,
         cellwidth=10,annotation_colors=annote.colors,
         filename = 'ciberSortTypes.pdf',width=18)
