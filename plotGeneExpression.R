##plot gene expression


require(reticulate)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
##get metadata file
#condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"


#reticulate::use_condaenv(condaenv)
syn=reticulate::import('synapseclient')
sync=syn$login()
return(sync)
##get annotations
tabres = sync$tableQuery('SELECT zScore,specimenID,Symbol,individualID,experimentalCondition FROM syn22878645')
tabres2 = sync$tableQuery('SELECT zScore,specimenID,Symbol,individualID,experimentalCondition FROM syn21222341')

rnaseq=rbind(tabres$asDataFrame(), tabres2$asDataFrame())%>%
  distinct()%>%
  mutate(specimenID=str_replace_all(specimenID,fixed('NF0007-2-M+C'),'NF0007-2- M+C'))%>%
  mutate(specimenID=str_replace_all(specimenID,fixed('NF0007-2-M+C+F'),'NF0007-2- M+C+F'))

annotes<-rnaseq%>%
  dplyr::select(specimenID,individualID,experimentalCondition)%>%
  separate(specimenID,into=c('altID','extra'),sep=' ',remove=F)%>%
  mutate(altID=stringr::str_replace(altID,'-$',''))%>%select(-extra)%>%
#  mutate(altID=stringr::str_replace(altID,'[-2-M+C]+','-2'))%>%
#  mutate(altID=stringr::str_replace(altID,'[007-2+C]+','007-2'))%>%
  distinct()

rownames(annotes)<-c()
annotes<-annotes%>%tibble::column_to_rownames('specimenID')


##update annote
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
nannotes<-apply(annotes,2,as.character)
rownames(nannotes)<-rownames(annotes)

mat<-rnaseq%>%
  dplyr::select(specimenID,zScore,Symbol)%>%
  tidyr::pivot_wider(values_from=zScore,names_from=specimenID,
                     values_fn=list(zScore=mean),values_fill=list(zScore=0))%>%
  tibble::column_to_rownames('Symbol')%>%as.matrix()

vars<-apply(mat,1,var,na.rm=T)%>%sort(decreasing=T)
  


#' Old plot using clusterProfiler
#' @export 
#' @import org.Hs.eg.db
#' @import clusterProfiler
plotOldGSEA<-function(genes.with.values,prefix,gsea_FDR=0.05){
  require(org.Hs.eg.db)
  # mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
  #    dplyr::rename(Gene='alias_symbol')
  
  #  genes.with.values<-genes.with.values%>%
  #    dplyr::left_join(mapping,by='Gene')%>%
  #    arrange(desc(value))
  genes.with.values<-arrange(genes.with.values,desc(value))
  # print(head(genes.with.values))
  genelist=genes.with.values$value
  names(genelist)=genes.with.values$Gene
  print(head(genelist))
  genelist<-sort(genelist,decreasing=TRUE)
  
  gr<-clusterProfiler::gseGO(unlist(genelist),ont="BP",keyType="SYMBOL",
                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = gsea_FDR)#,eps=1e-10)
  
  res<-filter(as.data.frame(gr),p.adjust<gsea_FDR)
  if(nrow(res)==0)
    return(gr)
  
  all_gseaGO<-res %>% 
    dplyr::rename(pathway = 'Description') %>% 
    arrange(NES) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up",
                                     NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    group_by(status) %>% 
    top_n(20, wt = abs(NES)) %>% 
    ungroup() %>% 
    ggplot2::ggplot(aes(x=reorder(pathway, NES), y=NES)) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
    coord_flip() +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") +
    labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
    ggtitle(paste('All',prefix))
  ggsave(paste0("allRegProts_", prefix,"_gseaGO_plot.pdf"), all_gseaGO, height = 8.5, width = 11, units = "in")
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("GO Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_GO.pdf'),width=10,height=10)
  
  df<-as.data.frame(gr)%>%mutate(Condition=prefix)
  return(df)
}


library(pheatmap)
#pheatmap(mat[names(vars)[1:1000],],annotation_col=annotes,cellwidth=10)

##which genes are differentially expressed between all 'treated' and 'none'?
#'
#'limmaTwoFactorDEAnalysis
#'uses Osama's code to compute de from limma
#'@author Osama
#'@import limma
#'@param data matrix
#'@param group1 ids
#'@param group2 ids
limmaTwoFactorDEAnalysis <- function(dat, sampleIDs.group1, sampleIDs.group2) {
  # Conduct DE expression analysis using limma from the expression matrix dat (group2 vs group1, group1 is reference)
  #
  # Args:
  #   dat: Expression data matrix, rows are genes, columns are samples
  #   sampleIDs.group1: Vector with ids of samples in reference group (eg. normal samples)
  #   sampleIDs.group2: Vector with ids of samples in interest group (eg. tumor samples) 
  #
  # Returns:
  #   limma Differential Expression results.
  #
  #http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html
  #http://genomicsclass.github.io/book/pages/using_limma.html
  #https://wiki.bits.vib.be/index.php/Tutorial:_Testing_for_differential_expression_I
  library(limma)
  fac <- factor(rep(c(2,1), c(length(sampleIDs.group2), length(sampleIDs.group1))))
  design <- model.matrix(~fac)
  fit <- lmFit(dat[,c(sampleIDs.group2, sampleIDs.group1)], design)
  fit <- eBayes(fit)
  print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="none")
  res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  return(arrange(res,P.Value))
}

plotDifferencesInCondition<-function(cond='None'){
  print(cond)
  cols=rownames(nannotes)[which(nannotes[,cond]=='TRUE')]
  ncols=setdiff(rownames(nannotes),cols)
  
  res=limmaTwoFactorDEAnalysis(mat,cols,ncols)%>%dplyr::arrange(P.Value)
  numGenes= nrow(subset(res,adj.P.Val<0.05))
  print(numGenes)
  print(nannotes)
  pheatmap(mat[res$featureID[1:50],],cellheight=10, cellwidth=10,
           annotation_col = as.data.frame(nannotes,stringsAsFactors = FALSE)%>%dplyr::select(-experimentalCondition),
           filename=paste0(cond,'top50transcripts.pdf'))
  ##we should probably do some sort of GSEA enrichment poltting as well.
  print(head(res))
  #nres<-res%>%dplyr::arrange(as.numeric(logFC))
  
  #genes<-nres$logFC
  #names(genes)<-nres$featureID
  res%>%dplyr::select(featureID,logFC)%>%dplyr::rename(value='logFC',Gene='featureID')%>%plotOldGSEA(.,cond)
  
}

library(ggfortify)
p1<-autoplot(prcomp(t(mat)),data=nannotes[colnames(mat),],shape='Media',colour='individualID')


ggsave('pcaOfAllSamples.png',p1)
#sapply(c('DMEM','Mammo','StemPro'),plotDifferencesInCondition)

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


orgs<-setdiff(annotes$altID,pats)

dlist<-lapply(orgs,function(pat){
  ##tannotes<-subset(annotes, individualID%in%c(pat,pats))%>%
    #dplyr::select(-individualID)
  
  iannote<-subset(annotes,altID==pat)%>%
    dplyr::select(experimentalCondition,Media,Cytokines,Forskolin='Forskoline')

  norm=iannote%>%subset(experimentalCondition=='None')%>%rownames()
  
  norcors<-sapply(setdiff(rownames(iannote),norm),function(x) cor(mat[,norm],mat[,x],method='spearman'))
 
  pdat<-iannote%>%subset(experimentalCondition!="None")%>%
    dplyr::select(Media,Cytokines,Forskolin)%>%
    cbind(Similarity=norcors)
  
  return(pdat)
})

names(dlist)<-orgs
ddf<-do.call(rbind,lapply(names(dlist),function(x) data.frame(Patient=x,dlist[[x]])))%>%
  tibble::rownames_to_column('altID')

write.csv(ddf,'orgCorrelations.csv',row.names=F)

names(dlist)<-orgs
plist<-lapply(orgs,function(pat){
  pdat<-dlist[[pat]]
  pdat%>%ggplot(aes(y=Similarity,x=Media,shape=Forskolin,color=Cytokines))+
    geom_point(aes(size=10))+
    ggtitle(pat)
})

library(cowplot)
res=cowplot::plot_grid(plotlist=plist)
ggsave(filename='corPlots.pdf',res,width=10)


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
  ggplot(aes(x=individualID,y=Similarity,fill=Organoid))+geom_boxplot()+scale_fill_viridis_d()
ggsave('cNFPatients.pdf',p2,width=10)


