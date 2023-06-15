##plot gene expression

source('loadOrganoidData.R')


library(ggfortify)
library(umap)
library(cowplot)

### Plot PCA of expression data
p1<-autoplot(prcomp(t(mat)),data=nannotes[colnames(mat),],shape='Media',colour='individualID')
#p2<-umap(t(mat), data=nannotes[colnames(mat),])

ggsave('pcaOfAllSamples.png',p1)

#sync$store(syn$File('pcaOfAllSamples.png',parentId='syn24827084'))


patindiv<-annotes[pats,'individualID']
specs<-rownames(annotes)


npannotes<-data.frame(individualID=c(annotes$individualID,gsub('tumor[0-9]*','',pats)),
             Media=c(annotes$Media,rep('Tumor',length(pats))),
             Cytokines=as.character(c(annotes$Cytokines,rep(FALSE,length(pats)))),
                           Forksoline=as.character(c(annotes$Forskoline,rep(FALSE,length(pats)))),
            cohort=c(ifelse(annotes$Media=='Tumor','cNF','Organoid'),rep('Biobank',length(pats))))

levels(npannotes$Media)<-levels(annotes$Media)
rownames(npannotes)<-c(rownames(annotes),pats)


annote.colors<-lapply(names(npannotes),function(x) c(`FALSE`='white',`TRUE`='darkgrey'))
names(annote.colors)<-names(npannotes)
annote.colors$Media<-media_pal
annote.colors$cohort<-c(cNF='white',Organoid='darkgrey',Biobank='grey')
annote.colors$individualID <- patient_pal

pheatmap(cor(mat,method='spearman'),
         annotation_col = npannotes,#%>%
           #dplyr::select(-individualID),
         annotation_row=npannotes,#%>%
         clustering_distance_rows = 'correlation',
         clustering_distance_cols='correlation',
          # dplyr::select(-individualID),
          cellheight=10,cellwidth = 10, annotation_colors=annote.colors,
         filename=paste0('heatmapOfAllCorrelations.pdf'))
#sync$store(syn$File('heatmapOfAllCorrelations.pdf',parentId='syn24827084'))


mat2<-rnaseq%>%
  subset(specimenID%in%specs)%>%
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  as.data.frame()%>%
  tibble::column_to_rownames('Symbol')%>%
  as.matrix()


shared<-intersect(colnames(mat),rownames(annotes))
ddf<-plotCorrelationBetweenSamps(mat[,shared],annotes[shared,],'geneExpression')

#sync$store(syn$File('geneExpressioncorPlots.pdf',parentId='syn24827084'))
write.table(ddf,file='tmp.csv',sep=',',row.names=F)
#sync$store(syn$build_table('RNAseq-based sample correlations','syn11374354','tmp.csv'))

badsamps<-subset(ddf,Similarity<0.6)$altID
print(paste('We have',length(badsamps),'bad samples we are removing'))

##now let's plot correlation with CNFs
nfsamps<-colnames(mat2)[grep('NF',colnames(mat))]
others<-setdiff(colnames(mat),nfsamps)
restab=do.call(rbind,lapply(nfsamps,function(x) cor(mat[,x],mat[,others],method='spearman')))
rownames(restab)<-nfsamps

cnfs<-annotes%>%
  dplyr::select(individualID,Media,extras)%>%
  as.data.frame()%>%
  tibble::rownames_to_column('cNF Sample')

p2<-restab%>%as.data.frame()%>%tibble::rownames_to_column('cNF Sample')%>%
  tidyr::pivot_longer(others,names_to='patient',values_to='Similarity')%>%
  left_join(cnfs)%>%
  mutate(`Biobank Patient`=stringr::str_replace_all(patient,'tumor[0-9]*',''))%>%
  mutate(`Biobank Patient`=stringr::str_replace_all(`Biobank Patient`,'patient','Patient '))%>%
  mutate(`Biobank Patient`=factor(`Biobank Patient`,levels=rev(c('Patient 1','Patient 2','Patient 3','Patient 4','Patient 5','Patient 6',
                                                                 'Patient 8','Patient 9','Patient 10','Patient 11','Patient 13'))))
p2%>%group_by(Media)%>%summarize(medSim=median(Similarity))

p3<-p2%>%
  rowwise()%>%
  mutate(extras=ifelse(Media=="Tumor"&&extras=="None","Tumor",extras))%>%
  ggplot(aes(x=`Biobank Patient`,y=Similarity,fill=Media))+geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=media_pal)+
  facet_grid(~extras)+coord_flip()

ggsave('cNFPatients.pdf',p2,width=10)

pat2=grep('NF0002',rownames(restab))

p2<-restab[-pat2,]%>%
  as.data.frame()%>%tibble::rownames_to_column('cNF Sample')%>%
  tidyr::pivot_longer(others,names_to='patient',values_to='Similarity')%>%
  left_join(cnfs)%>%
  mutate(`Biobank Patient`=stringr::str_replace_all(patient,'tumor[0-9]*',''))%>%
  mutate(`Biobank Patient`=stringr::str_replace_all(`Biobank Patient`,'patient','Patient '))%>%
  mutate(`Biobank Patient`=factor(`Biobank Patient`,levels=rev(c('Patient 1','Patient 2','Patient 3','Patient 4','Patient 5','Patient 6',
                                                                 'Patient 8','Patient 9','Patient 10','Patient 11','Patient 13'))))%>%
  rowwise()%>%
  mutate(extras=ifelse(Media=="Tumor"&&extras=="None","Tumor",extras))%>%
  ggplot(aes(x=`Biobank Patient`,y=Similarity,fill=Media))+geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=media_pal)+
  facet_grid(~extras)+coord_flip()

ggsave('cNFPatients_no0002.pdf',p2,width=10)


#ggsave('cNFPatients.png',p2,width=10)
##let's plot NF1 only?

nf1_counts<-rnaseq|>
  subset(Symbol=='NF1')|>
  left_join(annotes)|>
  subset(!is.na(experimentalCondition))|>
  ggplot(aes(x=individualID,fill=Media,shape=Cytokines,y=totalCounts))+
  geom_bar(stat='identity',position='dodge')+
  scale_fill_manual(values=media_pal)+
  facet_grid(Cytokines~.)

nf1_z<-rnaseq|>
  subset(Symbol=='NF1')|>
  left_join(annotes)|>
  subset(!is.na(experimentalCondition))|>
  ggplot(aes(x=individualID,fill=Media,shape=Cytokines,y=zScore))+
  geom_bar(stat='identity',position='dodge')+
  scale_fill_manual(values=media_pal)+
  facet_grid(Cytokines~.)

fp<-cowplot::plot_grid(nf1_counts,nf1_z)
ggsave('rnaseqcounts.pdf',fp)

pc<-tabres|>
  subset(experimentalCondition=='None')|>
  unique()|>
  left_join(annotes)|>
  mutate(isNF=(ifelse(Symbol=="NF1",TRUE,FALSE)))|>
  ggplot(aes(x=Symbol,y=totalCounts+0.1,col=isNF),xlabels=rep("",length(unique(rnaseq$Symbol))))+
  geom_point()+
  scale_y_log10()+
 # geom_bar(position='dodge',stat='identity')+
  facet_grid(specimenID~.)
  

