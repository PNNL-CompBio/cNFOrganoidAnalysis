##plottingFunctions


###here we standardized ploting/color schemes
###colors
if(!require(nationalparkcolors)){
  devtools::install_github("katiejolly/nationalparkcolors")
  library(nationalparkcolors)
}
library(pheatmap)
library(ggplot2)
if(!require(wesanderson)){
  install.packages("wesanderson")
  library(wesanderson)
}


##sort through color schemes to keep consistent throughout
pal = c(wesanderson::wes_palette("Chevalier1"),
        wesanderson::wes_palette("Cavalcanti1"),
        wesanderson::wes_palette("Zissou1"))

media_pal = wesanderson::wes_palette('Darjeeling1')[1:4]
names(media_pal)<-c('Mammo','DMEM','StemPro','Tumor')

org_pal = c("#212155","#2B3A8D","#3375B7","#6EC6EA","#CDE8F4")
names(org_pal)<-c("NF0009","NF0012","NF0002","NF0007","NF0008")

ctf_pal <- c("#27194D","#453286","#5D509D","#7165A9","#938ABE","#BFB7D9","#D1B4D3",
             "#BF8CBA","#AE5E9F","#873180","#572455")

names(ctf_pal)<-c("patient10","patient11","patient13","patient1","patient2","patient3","patient4","patient5","patient6","patient8","patient9")
patient_pal <- c(org_pal,ctf_pal)



#' getDifferencesInCOndition
#' Filters for a particular condition and data type
#' Computes limma differences at a p-value <0.05
#' @param mat
#' @param biga
#' @param cond
#' @param dataType
#' @param doPlot
#' @return number of genes
getDifferencesInCondition<-function(mat,biga,cond='None',dataType='geneExpression',doPlot=TRUE){
  print(cond)
  p1<-biga%>%dplyr::rename(val=cond)%>%
    subset(val==TRUE)
  p2<-biga%>%subset(None==TRUE)
  
  
  res=limmaTwoFactorDEAnalysis(mat,p1$specimenID,p2$specimenID)%>%
    dplyr::arrange(P.Value)
  mres<-subset(res,adj.P.Val<0.05)
  numGenes= nrow(mres)
  
  if(numGenes>1 && doPlot){
    pheatmap(mat[mres$featureID[1:min(50,numGenes)],biga$specimenID],cellheight=10, cellwidth=10,
             annotation_col = as.data.frame(nannotes,stringsAsFactors = FALSE)%>%dplyr::select(-experimentalCondition),
             filename=paste0(cond,'top50',dataType,'.pdf'))
  }
  ##we should probably do some sort of GSEA enrichment poltting as well.
  #print(head(res))
  #nres<-res%>%dplyr::arrange(as.numeric(logFC))
  
  return(numGenes)
}

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




#' explicitly plot correlations between samples from the same individual
plotCorrelationBetweenSamps<-function(mat,sannotes,prefix='geneExpression'){
  
  patids<-grep('patient',rownames(sannotes))
  
  samps<-sannotes%>%
    subset(Media%in%c('Tumor','None'))%>%rownames()|>
    setdiff(patids)
  
  #print(samps)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  ##now compute the correlation values
  dlist<-lapply(samps,function(norm){
  #  print(norm)
    iid<-sannotes[norm,'individualID']
    
    others<-setdiff(rownames(subset(sannotes,individualID==iid)),norm)
    others <-intersect(others,colnames(mat))
  #  norm <- intersect(rownames(sannotes)[grep(norm,rownames(sannotes))],colnames(mat))
   # print(others)
    norcors<-sapply(others,function(x) {
      cor(mat[,norm],mat[,x],method='spearman',use='pairwise.complete.obs')})
  #  print(norcors)
    pdat<-sannotes[others,]%>%#subset(Media!="None")%>%
      dplyr::select(Media,Cytokines,Forskolin)%>%
      cbind(Similarity=norcors)%>%
      tidyr::replace_na(list(Similarity=0.0))
   # print(pdat)
    return(pdat)
  })
  
  names(dlist)<-samps
  
  ddf<-do.call(rbind,lapply(names(dlist),function(x) data.frame(Patient=x,dlist[[x]])))%>%
    as.data.frame()%>%
    tibble::rownames_to_column('altID')%>%
 #   left_join(nann)%>%
    mutate(dataType=prefix)
  
  #write.csv(ddf,paste0(prefix,'orgCorrelations.csv'),row.names=F)

  names(dlist)<-samps
  plist<-lapply(samps,function(pat){
    pdat<-dlist[[pat]]
    pdat%>%ggplot(aes(y=Similarity,x=Media,shape=Forskolin,color=Cytokines))+
      geom_point(aes(size=10))+scale_colour_manual(values=pal)+
      ggtitle(pat)
  })
  
  #patCors<-lapply(samps,function(norm){
  #  print()
  #})
  print(head(ddf))
  res<-ggplot(ddf,aes(x=`Patient`,y=Similarity,col=Media))+geom_point()+
    scale_color_manual(values=pal)+
    facet_grid(Cytokines~Forskolin)
   # coord_flip()
  library(cowplot)
  ggsave(filename=paste0(prefix,'corPlots.pdf'),res,width=10)
  return(ddf)
  
  
}

#explicitly plot 
plotPatientCors<-function(mat,pannotes,nannotes,prefix='geneExpression'){
  patids<-grep('patient',pannotes$individualID)
  
  corTab<-cor(mat,method='spearman')%>%
    as.data.frame()%>%tibble::rownames_to_column('Sample')%>%
    pivot_longer(grep('patient',colnames(mat)),values_to='Correlation',names_to='Patient')%>%
    select(Sample,Patient,Correlation)
  
  patsonly<-grep('patient',corTab$Patient)
  corTab<-corTab[patsonly,]
  ##add in the cNF patient ids
  corTab<-corTab%>%
    left_join(tibble::rownames_to_column(pannotes,'Patient'))%>%
    select(Sample,Patient,Correlation,individualID)%>%distinct()%>%
    subset(Sample%in%orgs$specimenID)
  
  #now add in the media conditions
  corTab<-corTab%>%left_join(tibble::rownames_to_column(nannotes,'Sample'),by='Sample')%>%
    select(Sample,Patient,Correlation,`Biobank Patient`='individualID.x',Media,experimentalCondition)%>%
    distinct()%>%
    mutate(Additives=stringr::str_replace_all(experimentalCondition,"None|,*DMEM|,*Mammo|,*StemPro",""))%>%
    mutate(dataType=prefix)
  
  p<-ggplot(corTab,aes(x=`Biobank Patient`,y=Correlation,fill=Media))+geom_boxplot()+
    scale_fill_manual(values=pal)+facet_grid(~Additives)+coord_flip()
  ggsave(paste0(prefix,'patientCor.pdf'),p,width=10)
  return(corTab)
}