##run PLIER analysis on data

library(PLIER)
library(dplyr)
library(tidyr)
library(stringr)
#rc2path = "https://ndownloader.figshare.com/files/10881866"

##download these behemoth separately
#https://figshare.com/articles/recount_rpkm_RData/5716033/4

source("loadOrganoidData.R")
#then we want to looad this file
plier.results <- readRDS("../amlresistancenetworks/results/beatAMLpatientProfiling/recount_PLIER_model.RDS")

#'
#'Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD, Jaffe AE, Langmead B and Leek JT (2017). "Reproducible RNA-seq analysis using recount2." Nature Biotechnology. doi: 10.1038/nbt.3838
#'And the PLIER preprint, if you use the PLIER model:
#'  Mao W, Harmann B, Sealfon SC, Zaslavsky E, and Chikina M (2017). "Pathway-Level Information ExtractoR (PLIER) for gene expression data." bioRxiv. doi: 10.1101/116061
#'  

exprs.mat<-rnaseq%>%
  dplyr::select(specimenID,'Symbol',totalCounts)%>%
  pivot_wider(values_from=totalCounts, names_from=specimenID,
              values_fn=list(totalCounts=mean),
              values_fill=list(totalCounts=0.01))%>% #no zeroes allowed!
  tibble::column_to_rownames('Symbol')%>%as.matrix()

zvars<-which(apply(exprs.mat,1,var)==0)
if(length(zvars)>0)
  exprs.mat<-exprs.mat[-zvars,]

source("../multi-plier/util/plier_util.R")
pat.recount.b <- GetNewDataB(exprs.mat = exprs.mat,
                             plier.model = plier.results)

most.var<-apply(pat.recount.b,1,var,na.rm=T)
top.vars<-pat.recount.b[names(sort(most.var,decreasing=T))[1:20],]

###now that all the data are loaded we can compute the correlation, regression, and random forest values

lv.df<-pat.recount.b%>%
    as.data.frame()%>%
    tibble::rownames_to_column("Latent Variable")%>%
    pivot_longer(-`Latent Variable`,names_to='Sample',values_to="Loading")


lv.df<-dplyr::rename(lv.df,Gene='Latent Variable')

mat<-lv.df%>%
#  subset(specimenID%in%specs)%>%
  dplyr::select(Sample,Loading,Gene)%>%
  tidyr::pivot_wider(values_from=Loading,names_from=Sample,
                     values_fn=list(Loading=mean),values_fill=list(Loading=0))%>%
  tibble::column_to_rownames('Gene')%>%as.matrix()


##now compute how many genes are differentially expressed
condlist<-c('DMEM','StemPro','Cytokines','Mammo','Forskoline')
all.d<-lapply(condlist,function(x){
  plotDifferencesInCondition(mat,biga,x,'LVs')  
  
})
names(all.d)<-condlist
