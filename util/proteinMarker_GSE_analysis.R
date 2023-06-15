library(synapser)
library(synapserutils)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(tidyr)
library(stringr)

synLogin("","") #username,password

# queries two Synapse tables, conforms specimenID columns to standard
# table1 specimenID column data extracts relevant data from table2
query_tables <- function(table1synID,table2synID,query1cols,query2cols) {
  
  table1 = synTableQuery(toString(sprintf('SELECT %s FROM %s',query1cols,table1synID)))
  ihc_dat=as.data.frame(table1)
  ihc_dat$specimenID = gsub("- ","-",ihc_dat$specimenID)
  ihc_dat$specimenID = gsub(" ","-",ihc_dat$specimenID)
  table2 = synTableQuery(toString(sprintf('SELECT %s FROM %s',query2cols,table2synID)))
  rnaseq=as.data.frame(table2)
  rnaseq$specimenID = gsub("- ","-",rnaseq$specimenID)
  rnaseq$specimenID = gsub(" ","-",rnaseq$specimenID)
  
  combined<-left_join(ihc_dat,rnaseq,by=c("specimenID"))
  return(combined)
}

master_table = query_tables('syn23667394','syn22878645','*','Symbol,zScore,specimenID,individualID')


#Calculate correlation
calc_correlation <- function(subset_table,stat_method){
  cor.res<-subset_table%>%
    group_by(Symbol)%>% ## for each marker
    summarize(numSamps=n(),corVal=cor(Percent_DAB_Stained,zScore,method=stat_method))%>%
    arrange(desc(corVal))
  sigs<-subset(cor.res,abs(corVal)>0.65)
  return(sigs)
}

#Get table subset by proteinMarker
generate_ranked_genelist <- function(subset_table,stat_method,prot){
  sigs = calc_correlation(subset_table,stat_method)
  write.table(sigs,file=toString(sprintf('ranked_list_%s.txt',prot)))
  return(sigs)
}

library(org.Hs.eg.db)

do_gseGO_analysis <- function(sigs,prot){
  genelist = sigs$corVal
  names(genelist) = sigs$Symbol
  gse = gseGO(geneList=genelist,ont="ALL",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.05)
  write.table(gse,toString(sprintf("gsego_%s.tsv",prot)),sep="\t")
  return(gse)
}

gene_set_analysis <- function(master_table) {
  for (prot in unique(master_table$proteinMarker)) {
    subset_table = master_table[which(master_table$proteinMarker == prot),]
    sigs = generate_ranked_genelist(subset_table,'pearson',prot)
    gse = do_gseGO_analysis(sigs,prot)
    ## gse defined for further plotting analysis
  }
}

gene_set_analysis(master_table)

