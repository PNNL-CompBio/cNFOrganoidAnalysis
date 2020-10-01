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

rnaseq=rbind(tabres$asDataFrame(), tabres2$asDataFrame())

annotes<-rnaseq%>%
  dplyr::select(specimenID,individualID,experimentalCondition)%>%
  distinct()
rownames(annotes)<-c()
annotes<-annotes%>%tibble::column_to_rownames('specimenID')


##update annotes
for(x in c('kines','Mammo','DMEM','Forskoline','StemPro')){
  annotes[[x]]<-FALSE
  annotes[[x]][grep(x,annotes$experimentalCondition)]<-TRUE
 # annotes[[x]]=as.character(annotes)
}
annotes<-annotes%>%dplyr::rename(Cytokines='kines')#%>%dplyr::select(-experimentalCondition)
nannotes<-apply(annotes,2,as.character)
rownames(nannotes)<-rownames(annotes)

mat<-rnaseq%>%dplyr::select(specimenID,zScore,Symbol)%>%distinct()%>%
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
p1<-autoplot(prcomp(t(mat)),data=annotes,colour='experimentalCondition',shape='individualID')
pats=annotes$individualID[grep('patient',annotes$individualID)]


ggsave('pcaOfAllSamples.png',p1)
sapply(c('Cytokines','DMEM','Forskoline','Mammo','StemPro'),plotDifferencesInCondition)

sannotes<-nannotes%>%as.data.frame(stringsAsFactors=FALSE)%>%
  dplyr::select(DMEM,Mammo,Forskoline,Cytokines,StemPro,individualID)%>%
  dplyr::mutate(cohort=ifelse(individualID%in%pats,'cNF','Organoid'))
rownames(sannotes)<-rownames(nannotes)

pheatmap(cor(mat,method='spearman'),
         annotation_col = sannotes%>%
           dplyr::select(-individualID),
         annotation_row=sannotes%>%
           dplyr::select(-individualID),
          cellheight=10,cellwidth = 10,
         filename=paste0('heatmapOfAllConditions.pdf'))


for(pat in unique(annotes$individualID)){
  tannotes<-subset(sannotes, individualID%in%c(pat,pats))%>%
    dplyr::select(-individualID)
  #tannotes2<-tannotes%>%
  #  mutate(dataset=ifelse(individualID%in%pats,'cNF','Organoid'))%>%
  #  dplyr::select(-individualID)
  #rownames(tannotes2)<-rownames(tannotes)
      #dplyr::select(-c(individualID))
  pmat<-mat[,rownames(tannotes)]
  pmat
  pheatmap(cor(pmat,method='spearman'),
          annotation_col = tannotes,
          annotation_row=tannotes,
          cellheight=10,cellwidth = 10,
           filename=paste0('heatmapOf',pat,'Conditions.pdf'))
}


