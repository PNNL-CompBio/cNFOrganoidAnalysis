##plot gene expression

source('loadOrganoidData.R')


library(ggfortify)
library(umap)


### Plot PCA of expression data
p1<-autoplot(prcomp(t(mat)),data=nannotes[colnames(mat),],shape='Media',colour='individualID')
#p2<-umap(t(mat), data=nannotes[colnames(mat),])
ggsave('pcaOfAllSamples.png',p1)

sync$store(syn$File('pcaOfAllSamples.png',parentId='syn24827084'))


patindiv<-annotes[pats,'individualID']

###plot heatmap of correlations
sannotes<-nannotes%>%as.data.frame(stringsAsFactors=FALSE)%>%
  dplyr::mutate(cohort=ifelse(individualID%in%patindiv,'cNF','Organoid'))%>%
  dplyr::select(Media,Cytokines,Forskoline,individualID,cohort)

rownames(sannotes)<-rownames(nannotes)

annote.colors<-lapply(names(sannotes),function(x) c(`FALSE`='white',`TRUE`='black'))
names(annote.colors)<-names(sannotes)
annote.colors$Media<-c(None='white',StemPro='black',Mammo='darkgrey',DMEM='lightgrey')
annote.colors$cohort<-c(cNF='white',Organoid='black')


pheatmap(cor(mat,method='spearman'),
         annotation_col = sannotes%>%
           dplyr::select(-individualID),
         annotation_row=sannotes%>%
           dplyr::select(-individualID),
          cellheight=10,cellwidth = 10, annotation_colors=annote.colors,
         filename=paste0('heatmapOfAllCorrelations.pdf'))
sync$store(syn$File('heatmapOfAllCorrelations.pdf',parentId='syn24827084'))


mat2<-rnaseq%>%
  subset(specimenID%in%specs)%>%
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  tibble::column_to_rownames('Symbol')%>%as.matrix()


shared<-intersect(colnames(mat),rownames(annotes))
ddf<-plotCorrelationBetweenSamps(mat[,shared],annotes[shared,],'geneExpression')

sync$store(syn$File('geneExpressioncorPlots.pdf',parentId='syn24827084'))
write.table(ddf,file='tmp.csv',sep=',',row.names=F)
sync$store(syn$build_table('RNAseq-based sample correlations','syn11374354','tmp.csv'))

badsamps<-subset(ddf,Similarity<0.6)$altID
print(paste('We have',length(badsamps),'bad samples we are removing'))

##now let's plot correlation with CNFs
nfsamps<-colnames(mat2)[grep('NF',colnames(mat2))]
others<-setdiff(colnames(mat2),nfsamps)
restab=do.call(rbind,lapply(nfsamps,function(x) cor(mat2[,x],mat2[,others],method='spearman')))
rownames(restab)<-nfsamps

cnfs<-annotes%>%
  dplyr::select(individualID)%>%
  as.data.frame()%>%
  tibble::rownames_to_column('cNF Sample')

#p2<-restab%>%as.data.frame()%>%tibble::rownames_to_column(Organoid')%>%
#  tidyr::pivot_longer(others,names_to='cNF Sample',values_to='Similarity')%>%
#  left_join(cnfs)%>%
#  ggplot(aes(x=individualID,y=Similarity,fill=Organoid))+geom_boxplot()+scale_fill_manual(values=pal)
#ggsave('cNFPatients.pdf',p2,width=10)


condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
all.d<-lapply(condlist,function(x){
  
  res<-getDifferencesInCondition(mat,subset(biga,!specimenID%in%badsamps),x,'geneExpression',doPlot=FALSE)  
  
})
names(all.d)<-condlist
diffStats<-data.frame(all.d)%>%mutate(data='rnaSeq')



