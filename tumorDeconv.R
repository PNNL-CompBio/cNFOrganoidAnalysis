##run tumor deconvolution on organoids


library(immunedeconv)  
source("loadOrganoidData.R")

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
 # pheatmap(log2(mtab+0.01),annotation_col=select(df,-specimenID),
#      cellheight = 10,cellwidth=10,
#      file=paste0(prefix,'_',method,'Preds.pdf'), height=10,
#    labels_col=rep(" ",ncol(mtab)))
  
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

condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
mat<-mc%>%pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
    tibble::column_to_rownames('cell_type')
  
all.d<-lapply(condlist,function(x){
  plotDifferencesInCondition(mat,biga,x,'mcpTypes')  
})
names(all.d)<-condlist
plotCorrelationBetweenSamps(mat,annotes,'mcpCounter')

mat<-xc%>%pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
    tibble::column_to_rownames('cell_type')
  
all.x<-lapply(condlist,function(x){
    plotDifferencesInCondition(mat,biga,x,'xcellTypes')  
})
names(all.x)<-condlist
##now do a per-cell type correlation
plotCorrelationBetweenSamps(mat,annotes,'xcell')
