###synapse tables that were preprocessed were not including some of patient 7 information, so instead of trying to
###rerun i created a new gene expressionl oading script

library(readxl)
library(stringr)
##get gene name metadata
library(biomaRt)##load before reticulate is called

###Step 1: transcript to gene name mapping
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
filters<-c('ensembl_transcript_id_version')
attrib=c('ensembl_transcript_id_version','hgnc_symbol')
##first we get only those genes that are on nuclear chromosomes
my_chr <- c(1:22,'X','Y')
mapping<-biomaRt::getBM(human,attributes=attrib,filters='chromosome_name',values=my_chr)


##Step 2: get protein coding genes
options(timeout=1000)
path='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz'
if(!file.exists('gencode.v29.transcripts.txt')){
  res<-download.file(path,destfile = basename(path),method='curl')
  R.utils::gunzip(basename(path),overwrite=T)
  system(paste("grep protein_coding",gsub(".gz","",basename(path)),"|cut -d '|' -f 6 |uniq > gencode.v29.transcripts.txt"))
}
genes=read.table('gencode.v29.transcripts.txt')

##Step 3: get synapse files from custom-manifest
pat_manifest='rnaSeqProcessing/rnaseqpatientids.csv'
org_manifest='rnaSeqProcessing/rnaseqorgids.csv'

##get metadata file
reticulate::use_virtualenv('r-reticulate')
syn<-reticulate::import('synapseclient')
sync<-syn$login()


getCountsFromMetadata<-function(manifest,sync,mapping,genes){
  #now get the counts files themselves
  metadata<-read.csv(manifest,header=T)
  dres<-apply(metadata,1,function(x){
    p<-read.table(sync$get(x[[1]])$path,sep='\t', header=T)
    p$specimenID=x[[3]]
    p$individualID=x[[4]]
    p$experimentalCondition=x[[5]]
    return(p)
     
  })

  tres<-do.call(rbind,dres)|>
    dplyr::rename(ensembl_transcript_id_version='Name')#|>
  
  wgn <- tres|>dplyr::left_join(mapping)|>
    subset(hgnc_symbol!='')
  
  red.tab<-subset(wgn,hgnc_symbol%in%genes[,1])
  
  ##now we can sum reads and z score
  sum.tab<-red.tab|>
    group_by(specimenID,individualID,hgnc_symbol)|>
    mutate(totalCounts=sum(NumReads))|>
    dplyr::select(-c(Length,EffectiveLength,ensembl_transcript_id_version,TPM,NumReads))|>
    distinct()
  
  sum.tab<-sum.tab|>
    group_by(specimenID,individualID)|>
    mutate(zScore=(totalCounts-mean(totalCounts+0.001,na.rm=T))/sd(totalCounts,na.rm=T))|>
    dplyr::rename(Symbol='hgnc_symbol')
  
  return(sum.tab)
}

orgtab<-getCountsFromMetadata(org_manifest,sync,mapping,genes)

pattab<-getCountsFromMetadata(pat_manifest,sync,mapping,genes)|>
  ungroup()|>
  dplyr::select(-specimenID)|>
  dplyr::rename(specimenID='individualID')|>
  dplyr::mutate(individualID=stringr::str_replace(specimenID,'tumor[0-9]*',''))|>
  mutate(experimentalCondition='None')


##remove bad species from annnotations and RNA seq
badSpecs=c("NF0009-1 M","NF0009-1-M+C","NF0002-8-19-M+C+F","NF0002-8-19-D+C",
           "NF0002-8-19-D+C+F","NF0002-8-19-S+C")

##now fix the annotations!
rnaseq<-rbind(orgtab, pattab)%>%
  subset(!specimenID%in%badSpecs)|>
  distinct()%>%
  mutate(specimenID=str_replace_all(specimenID,'NF0007-4-D$','NF0007-4 D'))%>%
  mutate(specimenID=str_replace_all(specimenID,'NF0007-4-S$','NF0007-4 S'))%>%
  mutate(specimenID=str_replace_all(specimenID,'NF0007-4-M$','NF0007-4 M'))%>%
  mutate(experimentalCondition=str_replace_all(experimentalCondition,'Cytkines,Mammo','Cytokines,Mammo'))|>
  ungroup()

rnaseq$Media<-sapply(rnaseq$experimentalCondition,function(x){
  ifelse(length(grep('Mammo',x))>0,'Mammo',
         ifelse(length(grep('StemPro',x))>0,'StemPro',
                ifelse(length(grep('DMEM',x))>0,'DMEM','Tumor')))})

rnaseq$Cytokines<-sapply(rnaseq$experimentalCondition,function(x)
  ifelse(length(grep('Cytokines',x))>0,ifelse(length(grep("Forskoline",x))>0,'Cytokines,Forskolin','Cytokines'),'None'))

write.table(rnaseq,file='newGeneEx.csv',sep=',',row.names=F,col.names=T)


###lastly generate teh annotatiopns we need
##regular annotations
pannotes<-rnaseq|>
  dplyr::select(specimenID,individualID,experimentalCondition,Media,Cytokines)|>
  distinct()%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('specimenID')

##binary annotations
bannotes<-pannotes|>
  dplyr::select(-experimentalCondition)
bannotes$Cytokines=rep(FALSE,nrow(bannotes))
bannotes$Cytokines[grep('Cytokines',pannotes$experimentalCondition)]<-TRUE
bannotes$Forskolin=rep(FALSE,nrow(bannotes))
bannotes$Forskolin[grep('Forskolin',pannotes$experimentalCondition)]<-TRUE
bannotes$Cytokines=as.character(bannotes$Cytokines)
bannotes$Forskolin=as.character(bannotes$Forskolin)



##Create an annotationt able that mathces our color scheme
annote.colors<-lapply(names(bannotes),function(x) c(`FALSE`='white',`TRUE`='darkgrey'))
names(annote.colors)<-names(bannotes)
annote.colors$Media<-media_pal
#annote.colors$cohort<-c(cNF='white',Organoid='darkgrey',Biobank='grey')
annote.colors$individualID <- patient_pal



##keep track of organoid
orgs<-unique(orgtab$specimenID)
pats<-unique(pattab$specimenID)


##matrisome

gset <- read.table(sync$get('syn26199051')$path,header=TRUE,sep='\t')[-1,]
