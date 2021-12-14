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

if(FALSE){
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

}

full.tab<-sync$tableQuery('select * from syn23667404')$asDataFrame()%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,' ',''))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'-*M$',' M'))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'-*D$',' D'))%>%
  mutate(specimenID=stringr::str_replace_all(specimenID,'-*S$',' S'))
  
#%>%
#  mutate(specimenID=stringr::str_replace_all(specimenID,'[D|M|S]',''))

condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
plot.annotes<-annotes%>%
  select(-c(experimentalCondition,extras))%>%
      mutate(Cytokines=as.character(Cytokines))%>%
     mutate(Forskoline=as.character(Forskoline))
#  biga%>%
#    select(c(Cytokines,Forskoline,DMEM,Mammo,StemPro,specimenID))%>%

#  mutate(DMEM=as.character(DMEM))%>%
#  mutate(Mammo=as.character(Mammo))%>%
#  mutate(StemPro=as.character(StemPro))%>%
#  tibble::column_to_rownames('specimenID')

#ipal<-pal[1:length(unique(annotes$individualID))]
#names(ipal)<-unique(annotes$individualID)
ipal<-patient_pal
annote.colors<-list(Media=media_pal,#c(Mammo='darkgrey',Tumor='white',StemPro='black',DMEM='lightgrey'),
                    Cytokines=c(`FALSE`='white',`TRUE`='grey'),
                    Forskoline=c(`FALSE`='white',`TRUE`='grey'))
#names(annote.colors)<-names(plot.annotes)
annote.colors$individualID <- c(patient_pal,pal)[1:length(unique(pannotes$individualID))]
names(annote.colors$individualID)<-unique(pannotes$individualID)

##MCP Counter
mat<-subset(full.tab,method=='mcp_counter')%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
    tibble::column_to_rownames('cell_type')

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',method='ward.D2',
         cellwidth=10,annotation_colors=annote.colors,filename = 'all_mcpTypes.pdf',width=18)

mat<-subset(full.tab,method=='mcp_counter')%>%
  subset(!specimenID%in%pats)%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,method='ward.D2',
         #clustering_distance_rows = 'correlation',clustering_distance_cols='correlation',
         cellwidth=10,annotation_colors=annote.colors,filename = 'mcpTypes.pdf',width=18)
#all.d<-lapply(condlist,function(x){
#  getDifferencesInCondition(mat,biga,x,'mcpTypes',doPlot=FALSE)  
#})

#names(all.d)<-condlist
#cellCors<-data.frame(all.d)%>%mutate(dataType='mcpCounterPreds')
mcp.cor<-plotCorrelationBetweenSamps(mat,annotes,'mcpCounter')

##CIBERSORT

mat<-subset(full.tab,method=='cibersort')%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

#all.c<-lapply(condlist,function(x){
#    getDifferencesInCondition(mat,biga,x,'Cibersort',doPlot=FALSE)  
#})
#names(all.c)<-condlist
#cellCors<-data.frame(all.c)%>%mutate(dataType='cibersortPreds')%>%rbind(cellCors)
cs.cor<-plotCorrelationBetweenSamps(mat,annotes,'cibersort')

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',method='ward.D2',
         cellwidth=10,annotation_colors=annote.colors,
         filename = 'all_ciberSortTypes.pdf',width=18)


mat<-subset(full.tab,method=='cibersort')%>%
  subset(!specimenID%in%pats)%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,method='ward.D2',
        #clustering_distance_rows = 'correlation',clustering_distance_cols='correlation',
         cellwidth=10,annotation_colors=annote.colors,
         filename = 'ciberSortTypes.pdf',width=18)
#XCELL

mat<-subset(full.tab,method=='xcell')%>%
  pivot_wider(-c(method,experimentalCondition),values_from=score,names_from=specimenID)%>%
  tibble::column_to_rownames('cell_type')

#all.c<-lapply(condlist,function(x){
#    getDifferencesInCondition(mat,biga,x,'Cibersort',doPlot=FALSE)  
#})
#names(all.c)<-condlist
#cellCors<-data.frame(all.c)%>%mutate(dataType='cibersortPreds')%>%rbind(cellCors)
xc.cor<-plotCorrelationBetweenSamps(mat,annotes,'xcell')

pheatmap(log10(1+mat),annotation_col=plot.annotes,cellheight=10,#clustering_distance_rows = 'correlation',
         #clustering_distance_cols='correlation',
         cellwidth=10,annotation_colors=annote.colors,method='ward.D2',
         filename = 'all_xcellTypes.pdf',width=18)

rbind(cs.cor,mcp.cor,xc.cor)%>%
  ggplot(aes(x=Media,y=Similarity,fill=dataType,shape=Cytokines))+geom_boxplot()+
  geom_jitter()+scale_fill_manual(values=pal)
