library(shiny)
library(survival)
library(survMisc)

# nowotwory = c('COADREAD', 'GBM', 'GBMLGG', 'HNSC', 'KIPAN', 'KIRC', 'LGG', 'LUAD', 'LUSC',
#               'OV', 'STES', 'THCA', 'UCEC')
nowotwory <- list("GBMLGG", "BRCA", "KIPAN")

dane = NULL
for (nowotwor in nowotwory)
{
  dane1 = read.table(paste('pvalue/', nowotwor, '_pvalue.txt', sep=''), h=T)
  dane1$gen = rownames(dane1)
  dane1$nowotwor = nowotwor
  dane = rbind(dane, dane1)
}
geny<-unique(dane$gen)

for(nowotwor in nowotwory)
{
  load(paste('RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep='')) #Wczytanie pacjentów.
  load(paste('RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep='')) #Wczytanie mutacji.
  
  nowotwor.clinical_sub<-get(paste(nowotwor, '.clinical', sep=""))[,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]
  
  nowotwor.mutations_sub<-get(paste(nowotwor, '.mutations', sep=""))[!is.na(get(paste(nowotwor, '.mutations', sep=""))$Hugo_Symbol),c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]
  nowotwor.mutations_sub$Tumor_Sample_Barcode<-tolower(substr(nowotwor.mutations_sub$Tumor_Sample_Barcode, 1, 12))
  nowotwor.mutations_sub$Variant_Classification <-
    ifelse(is.na(as.numeric(nowotwor.mutations_sub$Variant_Classification)), as.character(nowotwor.mutations_sub$Variant_Classification), NA)
  nowotwor.mutations_sub$czy_jest<-1
  
  rm(list=paste(nowotwor, '.clinical', sep=''))
  rm(list=paste(nowotwor, '.mutations', sep=''))
  
  nowotwor.cli_mut<-nowotwor.clinical_sub
  nowotwor.cli_mut$status<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),0,1)
  nowotwor.cli_mut$time<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),nowotwor.cli_mut$patient.days_to_last_followup,nowotwor.cli_mut$patient.days_to_death)
  
  
  for(gen in geny){
    nowotwor.mutations_sub2<-nowotwor.mutations_sub[nowotwor.mutations_sub$Hugo_Symbol==gen,2:4]
    colnames(nowotwor.mutations_sub2) <- c('Tumor_Sample_Barcode', 'Variant_Classification', gen)
    nowotwor.cli_mut_temp <- merge(nowotwor.cli_mut, nowotwor.mutations_sub2, by.y='Tumor_Sample_Barcode', by.x='patient.bcr_patient_barcode', all.x=T, all.y=F)
    k <- ncol(nowotwor.cli_mut_temp)
    nowotwor.cli_mut_temp[,k]<-ifelse(is.na(nowotwor.cli_mut_temp[,k]),0,nowotwor.cli_mut_temp[,k])
    if(!dir.exists('Zbiory')){
      dir.create('Zbiory')}
    
    if(!dir.exists(paste('Zbiory/', nowotwor, sep=""))){
      dir.create(paste('Zbiory/', nowotwor, sep=""))}
    
    write.table(nowotwor.cli_mut_temp, file = paste('Zbiory/', nowotwor, '/', nowotwor, '_', gen, '_cli_mut.txt', sep=""))
  }
  
  
}

