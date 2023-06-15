##run tumor deconvolution on organoids


library(immunedeconv)  
library(dplyr)
library(tidyr)
library(GS)
library(pheatmap)
source("orgPlottingFunctions.R")
#source("rnaSeqProcessing/loadExpFromCounts.R")
#source("loadOrganoidData.R")
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
  
  df<-tab|>
    as.data.frame()|>
    dplyr::select(c(specimenID,experimentalCondition))%>%
      unique()#%>%
  #    rename(study='studyName')
  rownames(df)<-df$specimenID
  #save as heatmap with metadata
  mtab<-res%>%
    as.data.frame()|>
    dplyr::select(-cell_type)
  rownames(mtab)<-res$cell_type
  ##now tidy up data to table
  td<-tidyr::gather(res,key="specimenID",value="score",-cell_type )%>%
    left_join(df,by='specimenID')
  td$method=method
  return(td)
}

#if(FALSE){
##NEED TO DOWNLOD FILES MANUALLY I CANNOT DISTRIBUTE
set_cibersort_binary('./CIBERSORT.R')
set_cibersort_mat('./LM22.txt')

xc<-runImmuneDeconv(rnaseq,'xcell')
mc<-runImmuneDeconv(rnaseq,'mcp_counter')
cs<-runImmuneDeconv(rnaseq,'epic')
#cs<-runImmuneDeconv(rnaseq,'cibersort')


full.tab<-rbind(xc,mc,cs)

write.table(full.tab,'deconvfile.tsv',sep='\t',quote=FALSE,row.names=F)
#tab<-syn$build_table('cNF Organoid Deconvolution','syn11374354','file.tsv')
#sync$store(tab)


##MCP Counter
mat<-subset(full.tab,method=='mcp_counter')%>%
  pivot_wider(id_cols=-c(method,experimentalCondition),
              values_from=score,names_from=specimenID)%>%
    tibble::column_to_rownames('cell_type')

pheatmap(log10(1+mat),annotation_col=bannotes,cellheight=10,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',method='ward.D2',
         cellwidth=10,annotation_colors=annote.colors,filename = 'all_mcpTypes.pdf',width=18)

mat<-subset(full.tab,method=='mcp_counter')%>%
  subset(!specimenID%in%pats)%>%
  pivot_wider(id_cols=-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

pheatmap(log10(1+mat),annotation_col=bannotes,cellheight=10,method='ward.D2',
         #clustering_distance_rows = 'correlation',clustering_distance_cols='correlation',
         cellwidth=10,annotation_colors=annote.colors,filename = 'mcpTypes.pdf',width=18)

#mcp.cor<-plotCorrelationBetweenSamps(mat,annotes,'mcpCounter')

##CIBERSORT

mat<-subset(full.tab,method=='epic')%>%
  pivot_wider(id_cols=-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

#all.c<-lapply(condlist,function(x){
#    getDifferencesInCondition(mat,biga,x,'Cibersort',doPlot=FALSE)  
#})
#names(all.c)<-condlist
#cellCors<-data.frame(all.c)%>%mutate(dataType='cibersortPreds')%>%rbind(cellCors)
#cs.cor<-plotCorrelationBetweenSamps(mat,annotes,'cibersort')

pheatmap(log10(1+mat),annotation_col=bannotes,cellheight=10,clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',method='ward.D2',
         cellwidth=10,annotation_colors=annote.colors,
         filename = 'all_epicTypes.pdf',width=18)


mat<-subset(full.tab,method=='epic')%>%
  subset(!specimenID%in%pats)%>%
  pivot_wider(id_cols=-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

pheatmap(log10(1+mat),annotation_col=bannotes,cellheight=10,method='ward.D2',
        #clustering_distance_rows = 'correlation',clustering_distance_cols='correlation',
         cellwidth=10,annotation_colors=annote.colors,
         filename = 'epicTypes.pdf',width=18)
#XCELL

mat<-subset(full.tab,method=='xcell')%>%
  subset(!specimenID%in%pats)%>%
  pivot_wider(id_cols=-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')



pheatmap(log10(1+mat),annotation_col=bannotes,cellheight=10,#clustering_distance_rows = 'correlation',
         #clustering_distance_cols='correlation',
         cellwidth=10,annotation_colors=annote.colors,method='ward.D2',
         filename = 'xcellTypes.pdf',width=18)


#rbind(cs.cor,mcp.cor,xc.cor)%>%
#  ggplot(aes(x=Media,y=Similarity,fill=dataType,shape=Cytokines))+geom_boxplot()+
#  geom_jitter()+scale_fill_manual(values=pal)

######last part of figure is matrisome

mat.rest <-rnaseq|>
  subset(!specimenID%in%pats)|>
  dplyr::select(specimenID,Symbol,totalCounts)|>
  dplyr::distinct()|>
  tidyr::pivot_wider(names_from='specimenID',values_from='totalCounts')|>
  tibble::column_to_rownames('Symbol')|>
  as.matrix()

GSET <- apply(gset,2,function(x) setdiff(x,""))
res <- gsva(mat.rest,GSET,method="ssgsea",min.sz=1,max.sz=Inf,mx.diff=TRUE,parallel.sz=1,verbose=TRUE)
#res <- t(as.data.frame(res))
print(res)


pheatmap(res,annotation_col = bannotes,annotation_colors=annote.colors,
         cellheight=10,cellwidth=10,filename='matrisomeGSVA.pdf')

