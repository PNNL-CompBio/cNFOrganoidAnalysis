##plot gene expression

#this file loads the ploting functions
source("orgPlottingFunctions.R")
source('loadExpFromCounts.R')


library(ggfortify)
library(umap)
library(cowplot)
library(dplyr)

###prelim_fig plot pca
mat <- rnaseq%>%
  dplyr::ungroup()|>
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  tibble::column_to_rownames('Symbol')%>%as.matrix()

### Plot PCA of expression data
p1<-autoplot(prcomp(t(mat)),data=pannotes[colnames(mat),],shape='Media',colour='individualID')
#p2<-umap(t(mat), data=nannotes[colnames(mat),])

ggsave('pcaOfAllSamples.png',p1)


##first let's remove patient 2
pat2=setdiff(rownames(bannotes),rownames(bannotes)[grep('NF0002',rownames(bannotes))])
bannotes<-bannotes[pat2,]
mat<-mat[,pat2]


pheatmap(cor(mat,method='pearson'),
         annotation_col = bannotes,#%>%
         annotation_row= bannotes,#%>%
         clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',
          # dplyr::select(-individualID),
          cellheight=10,cellwidth = 10, annotation_colors=annote.colors,
         filename=paste0('heatmapOfAllCorrelations.pdf'))


###now we can do just an organoid based clustering
omat<-rnaseq%>%
  subset(specimenID%in%intersect(pat2,orgs))%>%
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  as.data.frame()%>%
  tibble::column_to_rownames('Symbol')%>%
  as.matrix()

##create new annotations with just organoids
obannotes<-bannotes[-grep('patient',bannotes$individualID),]
shared<-intersect(colnames(omat),rownames(obannotes))
ddf<-plotCorrelationBetweenSamps(omat[,shared],obannotes[shared,],'Fig4A_')

#sync$store(syn$File('geneExpressioncorPlots.pdf',parentId='syn24827084'))
write.table(ddf,file='fig4a_data.csv',sep=',',row.names=F)
#sync$store(syn$build_table('counts-based sample correlations','syn11374354','tmp.csv'))

badsamps<-subset(ddf,Similarity<0.6)$altID
print(paste('We have',length(badsamps),'bad samples we are removing'))

##now let's plot correlation with CNFs
nfsamps<-colnames(mat)[grep('NF',colnames(mat))]
others<-setdiff(colnames(mat),nfsamps)
restab=do.call(rbind,lapply(nfsamps,function(x) cor(mat[,x],mat[,others],method='spearman')))
rownames(restab)<-nfsamps

cnfs<-bannotes%>%
  dplyr::select(individualID,Media,Cytokines,Forskolin)%>%
  as.data.frame()%>%
  tibble::rownames_to_column('cNF Sample')

p2<-restab%>%as.data.frame()%>%tibble::rownames_to_column('cNF Sample')%>%
  tidyr::pivot_longer(others,names_to='patient',values_to='Similarity')%>%
  left_join(cnfs)%>%
  dplyr::mutate(`Biobank Patient`=stringr::str_replace_all(patient,'tumor[0-9]*',''))%>%
  dplyr::mutate(`Biobank Patient`=stringr::str_replace_all(`Biobank Patient`,'patient','Patient '))%>%
  dplyr::mutate(`Biobank Patient`=factor(`Biobank Patient`,levels=c('Patient 1','Patient 2','Patient 3','Patient 4','Patient 5','Patient 6',
                                                                 'Patient 8','Patient 9','Patient 10','Patient 11','Patient 13')))
p2$Extras=rep("None",nrow(p2))
p2$Extras[which(p2$Media=='Tumor')]<-'Tumor'
p2$Extras[intersect(which(p2$Cytokines=='TRUE'),which(p2$Forskolin=='TRUE'))]<-'F+C'
p2$Extras[intersect(which(p2$Cytokines=='TRUE'),which(p2$Forskolin=='FALSE'))]<-'C only'
p2$Extras[intersect(which(p2$Cytokines=='FALSE'),which(p2$Forskolin=='TRUE'))]<-'F only'

p3<-p2%>%
  ggplot(aes(x=`Biobank Patient`,y=Similarity,fill=Media))+geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=media_pal)+
  facet_grid(Extras~.)+
  theme_classic()



p4<-p2%>%
  ggplot(aes(x=`Biobank Patient`,y=Similarity,fill=Media))+geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=media_pal)+
  facet_grid(Extras~.)+
  scale_y_continuous(limits=c(0.5,1))+
  theme_classic()

ggsave('fig4b_organoidCorrelation.pdf',p3,height=12,width=10)
ggsave('fig4b_organoidCorrelation_largerY.pdf',p4,height=12,width=10)

write.csv(p2,file='fig4b_data.csv',row.names=F)

##now let's do the deconvolution
