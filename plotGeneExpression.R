##plot gene expression

source('loadOrganoidData.R')


library(ggfortify)
p1<-autoplot(prcomp(t(mat)),data=nannotes[colnames(mat),],shape='Media',colour='individualID')
ggsave('pcaOfAllSamples.png',p1)


sannotes<-nannotes%>%as.data.frame(stringsAsFactors=FALSE)%>%
  dplyr::select(Media,Cytokines,Forskoline,individualID)%>%
  dplyr::mutate(cohort=ifelse(individualID%in%pats,'cNF','Organoid'))
rownames(sannotes)<-rownames(nannotes)

pheatmap(cor(mat,method='spearman'),
         annotation_col = sannotes%>%
           dplyr::select(-individualID),
         annotation_row=sannotes%>%
           dplyr::select(-individualID),
          cellheight=10,cellwidth = 10,
         filename=paste0('heatmapOfAllConditions.pdf'))



mat2<-rnaseq%>%
  subset(specimenID%in%specs)%>%
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  tibble::column_to_rownames('Symbol')%>%as.matrix()


##now let's plot correlation with CNFs
nfsamps<-colnames(mat2)[grep('NF',colnames(mat2))]
others<-setdiff(colnames(mat2),nfsamps)
restab=do.call(rbind,lapply(nfsamps,function(x) cor(mat2[,x],mat2[,others],method='spearman')))
rownames(restab)<-nfsamps

cnfs<-annotes%>%
  dplyr::select(individualID)%>%
  as.data.frame()%>%
  tibble::rownames_to_column('cNF Sample')

p2<-restab%>%as.data.frame()%>%tibble::rownames_to_column('Organoid')%>%
  tidyr::pivot_longer(others,names_to='cNF Sample',values_to='Similarity')%>%
  left_join(cnfs)%>%
  ggplot(aes(x=individualID,y=Similarity,fill=Organoid))+geom_boxplot()+scale_fill_manual(values=pal)
ggsave('cNFPatients.pdf',p2,width=10)


condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
all.d<-lapply(condlist,function(x){
  plotDifferencesInCondition(mat,biga,x,'geneExpression')  
  
})
names(all.d)<-condlist

plotCorrelationBetweenSamps(mat,annotes,'geneExpression')

