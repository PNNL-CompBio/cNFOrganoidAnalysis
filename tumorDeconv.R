##run tumor deconvolution on organoids


library(immunedeconv)  
source("loadOrganoidData.R")


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

condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
mat<-mc%>%pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
    tibble::column_to_rownames('cell_type')
  
plot.annotes<-biga%>%
    select(c(Cytokines,Forskoline,DMEM,Mammo,StemPro,specimenID))%>%
    mutate(Cytokines=as.character(Cytokines))%>%
  mutate(Forskoline=as.character(Forskoline))%>%
  mutate(DMEM=as.character(DMEM))%>%mutate(Mammo=as.character(Mammo))%>%mutate(StemPro=as.character(StemPro))%>%
  tibble::column_to_rownames('specimenID')

pheatmap(log10(1+mat),annotation_col=plot.annotes,filename = 'mcpTypes.pdf',width=10)

all.d<-lapply(condlist,function(x){
  getDifferencesInCondition(mat,biga,x,'mcpTypes',doPlots=FALSE)  
})
names(all.d)<-condlist
cellCors<-data.frame(all.d)%>%mutate(dataType='mcpCounterPreds')

ddf<-plotCorrelationBetweenSamps(mat,annotes,'mcpCounter')


mat<-cs%>%pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

all.c<-lapply(condlist,function(x){
    getDifferencesInCondition(mat,biga,x,'mcpCounterTypes',doPlots=FALSE)  
})

cellCors<-data.frame(all.c)%>%mutate(dataType='cibersortPreds')%>%rbind(cellCors)


pheatmap(log10(1+mat),annotation_col=plot.annotes,filename = 'ciberSortTypes.pdf',width=10)


all.x<-lapply(condlist,function(x){
    getDifferencesInCondition(mat,biga,x,'xcellTypes',doPlots=FALSE)  
})
names(all.x)<-condlist

cellCors<-data.frame(all.x)%>%mutate(dataType='xCellPreds')%>%
  rbind(cellCors)
  
##now do a per-cell type correlation
dcors<-plotCorrelationBetweenSamps(mat,annotes,'xcell')
