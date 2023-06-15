##plot gene expression

source('loadOrganoidData.R')


library(ggfortify)
library(umap)
library(readxl)

drug<-read_xlsx('drug_screening/nf2_7_4_9_12_heatmap_data_W-names.xlsx')|>
  tidyr::pivot_longer(cols=c(3,4,5,6),names_to='sample',values_to='viability')|>
  tidyr::separate(sample,into=c('sarc','num1','num2','num3','num4','num5'))|>
  dplyr::mutate(individualID=stringr::str_to_upper(stringr::str_replace(sarc,'sarc00','')))


drugvals<-drug|>
  dplyr::select(Drug='figure_name',individualID,viability)

drugvars<-drugvals|>
  group_by(Drug)|>
  dplyr::mutate(var=var(viability))|>
  dplyr::select(-c(individualID,viability))|>
  distinct()

sorted<-arrange(drugvars,desc(var))
print(sorted)

highvar<-subset(sorted,var>200)
print(paste('found',nrow(highvar),'drugs with variance over 200'))

library(ggridges)

drugPlots<-subset(drugvals,Drug%in%highvar$Drug)%>%
  ggplot(.,aes(y=viability,x=Drug,fill=individualID))+geom_bar(stat='identity',position = 'dodge')+
  scale_fill_viridis_d()

drugPlots

ggsave(drugPlots,file='variableDrugs.png')

###now that we have some variable drugs, we can begin to do differential expression
zscore<-rnaseq|>
  dplyr::select(zScore,specimenID,Symbol)|>
  tidyr::pivot_wider(names_from=specimenID,values_from='zScore')|>
  tibble::column_to_rownames('Symbol')

#zscore=do.call(cbind,apply(zscore,2,unlist))

tcounts<-rnaseq|>
  dplyr::select(totalCounts,specimenID,Symbol)|>
  tidyr::pivot_wider(names_from=specimenID,values_from='totalCounts')|>
  tibble::column_to_rownames('Symbol')

#tcounts=do.call(cbind,apply(tcounts,2,unlist))

dres<-lapply(unique(highvar$Drug),function(d){
    print(d)
    dvals<-subset(drugvals,Drug==d)%>%arrange(viability)
    botvals=dvals$individualID[c(1,2)]
    topvals=dvals$individualID[c(3,4)]
    #print(topvals)
    #print(botvals)
    topcols<-colnames(zscore)[c(grep(topvals[1],colnames(zscore)),grep(topvals[2],colnames(zscore)))]
    botcols<-colnames(zscore)[c(grep(botvals[1],colnames(zscore)),grep(botvals[2],colnames(zscore)))]
    
    res<-limmaTwoFactorDEAnalysis(as.matrix(zscore),topcols,botcols)
    res
})

names(dres)<-unique(highvar$Drug)

dcounts<-lapply(dres,function(x) x|>subset(adj.P.Val<0.05)|>subset(abs(logFC)>1)|>nrow())
counttab<-cbind(Drug=unique(highvar$Drug),dcounts)
print(counttab)


##plot drugs

allplots<-lapply(rownames(counttab),function(dname){
  genes<-dres[[dname]]|>
    subset(adj.P.Val<0.05)|>
    subset(abs(logFC)>1)|>
    dplyr::select(featureID)


  geneMat<-tcounts[genes$featureID,-grep('patient',colnames(tcounts))]
  
  annotes<-rnaseq|>dplyr::select(specimenID,individualID,experimentalCondition)|>distinct()|>
    left_join(drug)|>
    subset(figure_name==dname)|>
    dplyr::select(individualID,specimenID,viability)|>
    distinct()|>
    subset(!is.na(individualID))|>
    tibble::column_to_rownames('specimenID')

  library(pheatmap)
  pheatmap(log10(0.01+geneMat)[,rownames(annotes)],annotation_col=annotes,filename=paste0(dname,'diffexGenes.pdf'))
})
##obviously there are only so many ways we can sort the patients, but for now we can use the leapR library to do the enrichment
library(leapR)
data('krbpaths')

##format diffex in columns
logvals<-do.call(cbind,lapply(names(dres),function(x){
  mat=dres[[x]]|>
    dplyr::select(logFC)
  colnames(mat)[1]<-x
  mat
}))

pvals<-do.call(cbind,lapply(names(dres),function(x){
  mat=dres[[x]]|>
    dplyr::select(adj.P.Val)
  colnames(mat)[1]<-x
  mat
}))

##now we can do the analysis

res<-do.call(cbind,lapply(colnames(pvals),function(p){
  lr<-leapR::leapR(geneset=krbpaths,enrichment_method='enrichment_in_sets',pvals,
                   primary_columns=p,threshold=0.05,
                   greaterthan=FALSE)|>
    dplyr::select(BH_pvalue)
  colnames(lr)<-p
  lr
}))

nas<-which(apply(res,1,function(x) all(x==1)))
pheatmap(res[-nas,])

ores<-do.call(cbind,lapply(colnames(pvals),function(p){
  lr<-leapR::leapR(geneset=krbpaths,enrichment_method='enrichment_in_order',logvals,
                   primary_columns=p)|>
    dplyr::select(pvalue)
  colnames(lr)<-p
  lr
}))

nas<-which(apply(ores,1,function(x) all(is.na(x))))
pheatmap(ores[-nas,])
