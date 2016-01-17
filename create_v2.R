library(shiny)
library(survival)
library(survMisc)

nowotwory = c('STES')
#, 'COADREAD', 'GBM', 'GBMLGG', 'HNSC', 'KIPAN', 'KIRC', 'LGG', 'LUAD', 'LUSC', 'OV', 'STES', 'THCA', 'UCEC')

geny <- read.table("p_value/lista_interesujacych_genow.txt", h=T)

for(nowotwor in nowotwory){
  #definiujemy ścieżkę do zbiorów
  dir_zbiory <- 'C:/Users/Sebastian/OneDrive/Studia/II stopien/Warsztaty badawcze - pbiecek/'
  
  load(paste(dir_zbiory, 'RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep='')) #Wczytanie pacjentów.
  load(paste(dir_zbiory, 'RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep='')) #Wczytanie mutacji.
  
  nowotwor.clinical_sub<-get(paste(nowotwor, '.clinical', sep=""))[,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]
  
  nowotwor.mutations_sub<-get(paste(nowotwor, '.mutations', sep=""))[!is.na(get(paste(nowotwor, '.mutations', sep=""))$Hugo_Symbol),c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]
  nowotwor.mutations_sub$Tumor_Sample_Barcode<-tolower(substr(nowotwor.mutations_sub$Tumor_Sample_Barcode, 1, 12))
  nowotwor.mutations_sub$Variant_Classification <-
    ifelse(is.na(as.numeric(nowotwor.mutations_sub$Variant_Classification)), as.character(nowotwor.mutations_sub$Variant_Classification), NA)
  nowotwor.mutations_sub$czy_jest<-1
  
  rm(list=paste(nowotwor, '.clinical', sep=''))
  rm(list=paste(nowotwor, '.mutations', sep=''))
  
  ktory <- 1
  for(g in 1:nrow(geny)){
    gen <- geny[g,]
    print(paste(nowotwor, ':', ktory, '/', nrow(geny)))
    ktory <- ktory + 1
    
    nowotwor.cli_mut.I<-nowotwor.clinical_sub
    nowotwor.cli_mut.I$status<-ifelse(is.na(nowotwor.cli_mut.I$patient.days_to_death),0,1)
    nowotwor.cli_mut.I$time<-ifelse(is.na(nowotwor.cli_mut.I$patient.days_to_death),nowotwor.cli_mut.I$patient.days_to_last_followup,nowotwor.cli_mut.I$patient.days_to_death)
    
    nowotwor.cli_mut.V<-nowotwor.clinical_sub
    nowotwor.cli_mut.V$status<-ifelse(is.na(nowotwor.cli_mut.V$patient.days_to_death),0,1)
    nowotwor.cli_mut.V$time<-ifelse(is.na(nowotwor.cli_mut.V$patient.days_to_death),nowotwor.cli_mut.V$patient.days_to_last_followup,nowotwor.cli_mut.V$patient.days_to_death)
    
    nowotwor.mutations_sub_gen.I<-nowotwor.mutations_sub[nowotwor.mutations_sub$Hugo_Symbol==gen,c("Tumor_Sample_Barcode","czy_jest")]
    colnames(nowotwor.mutations_sub_gen.I) <- c('Tumor_Sample_Barcode', paste('I.',gen, sep=""))
    nowotwor.cli_mut.I <- merge(nowotwor.cli_mut.I, nowotwor.mutations_sub_gen.I, by.y='Tumor_Sample_Barcode', by.x='patient.bcr_patient_barcode', all.x=T, all.y=F)
    k <- ncol(nowotwor.cli_mut.I)
    nowotwor.cli_mut.I[,k]<-ifelse(is.na(nowotwor.cli_mut.I[,k]),0,nowotwor.cli_mut.I[,k])
    nowotwor.cli_mut.I <- unique(nowotwor.cli_mut.I)
    
    nowotwor.mutations_sub_gen.V<-nowotwor.mutations_sub[nowotwor.mutations_sub$Hugo_Symbol==gen,c("Tumor_Sample_Barcode","Variant_Classification")]
    colnames(nowotwor.mutations_sub_gen.V) <- c('Tumor_Sample_Barcode', paste('V.',gen, sep=""))
    nowotwor.cli_mut.V <- merge(nowotwor.cli_mut.V, nowotwor.mutations_sub_gen.V, by.y='Tumor_Sample_Barcode', by.x='patient.bcr_patient_barcode', all.x=T, all.y=F)
    k <- ncol(nowotwor.cli_mut.V)
    nowotwor.cli_mut.V[,k]<-ifelse(is.na(nowotwor.cli_mut.V[,k]),0,nowotwor.cli_mut.V[,k])
    
    if(!dir.exists('Zbiory_v2')){
      dir.create('Zbiory_v2')}
    
    if(!dir.exists(paste('Zbiory_v2/', nowotwor, sep=""))){
      dir.create(paste('Zbiory_v2/', nowotwor, sep=""))}
    
    write.table(nowotwor.cli_mut.I, file = paste('Zbiory_v2/', nowotwor, '/', nowotwor, '_', gen, '_I.txt', sep=""))
    write.table(nowotwor.cli_mut.V, file = paste('Zbiory_v2/', nowotwor, '/', nowotwor, '_', gen, '_V.txt', sep=""))
    
  }
  
}