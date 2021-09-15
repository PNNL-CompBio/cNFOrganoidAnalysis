Sys.setenv(RETICULATE_PYTHON = '/Users/bade228/opt/anaconda3/envs/r2/bin/python3')
library(reticulate)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(tidyr)
library(stringr)
library(pathview)
library(org.Hs.eg.db)
syn=reticulate::import('synapseclient')
sync=syn$login()

# queries two Synapse tables, conforms specimenID columns to standard
# table1 specimenID column data extracts relevant data from table2
query_tables <- function(table1synID,table2synID,query1cols,query2cols) {
  
  ihc_dat = sync$tableQuery(paste0('SELECT ',query1cols,' FROM ',table1synID))$asDataFrame()
  ihc_dat$specimenID = gsub(" ","-",ihc_dat$specimenID)
  rnaseq = sync$tableQuery(paste0('SELECT ',query2cols,' FROM ',table2synID))$asDataFrame()
  rnaseq$specimenID = gsub(" ","-",rnaseq$specimenID)
  combined<-left_join(ihc_dat,rnaseq,by=c("specimenID"))
  return(combined)
}

master_table = query_tables('syn24175711','syn22878645','specimenID,proteinMarker,Fraction_DAB_stained','Symbol,zScore,specimenID,individualID')
#write.table(master_table,'master_list.tsv',sep='\t')
#master_table <- read.table('master_list.tsv')

# +
#Calculate correlation
calc_correlation <- function(subset_table,stat_method){
  cor.res<-subset_table%>%
    group_by(Symbol)%>% ## for each marker
    summarize(numSamps=n(),corVal=cor(Fraction_DAB_stained,zScore,method=stat_method))%>%
    arrange(desc(corVal))
  sigs<-subset(cor.res,abs(corVal)>0.65)
  #sigs <- cor.res
  return(sigs)
}

#Get table subset by proteinMarker
generate_ranked_genelist <- function(subset_table,stat_method,prot){
  sigs = calc_correlation(subset_table,stat_method)
  #write.table(sigs,file=toString(sprintf('ranked_list_%s.txt',prot)))
  return(sigs)
}

do_gseGO_analysis <- function(sigs,prot){
  genelist = sigs$corVal
  names(genelist) = sigs$Symbol
  gse = gseGO(geneList=genelist,ont="ALL",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.05)
  #write.table(gse,toString(sprintf("gsego_%s.tsv",prot)),sep="\t")
  return(gse)
}

gene_set_analysis <- function(master_table) {
  #library(DOSE)
  for (prot in unique(master_table$proteinMarker)) {
    subset_table = master_table[which(master_table$proteinMarker == prot),]
    sigs = generate_ranked_genelist(subset_table,'pearson',prot)
    gse = do_gseGO_analysis(sigs,prot)
    print(gse)
    ## gse defined for further plotting analysis
    #dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
    #emapplot(gse,showCategory=10)
  }
}
# -

gene_set_analysis(master_table)

