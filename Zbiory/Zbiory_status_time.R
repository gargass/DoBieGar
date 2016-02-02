nowotwory <- c('BRCA', 'COADREAD', 'GBM', 'GBMLGG', 'HNSC', 'KIPAN', 'KIRC', 'LGG', 'LUAD', 'LUSC', 'OV', 'STES', 'THCA', 'UCEC')


dir_zbiory <- 'C:/Users/Gosia/Desktop/Biecek/'


for(nowotwor in nowotwory){
  print(nowotwor)
load(paste(dir_zbiory, 'RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep='')) #Wczytanie pacjentów.
#load(paste(dir_zbiory, 'RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep='')) #Wczytanie mutacji.

nowotwor.geny <- read.table(paste('Zbiory/', nowotwor, '.txt', sep=""))
nowotwor.geny <- nowotwor.geny[
  order(nowotwor.geny$pacjent), ]

nowotwor.clinical <-get(paste(nowotwor, '.clinical', sep=""))[
  ,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]

rm(list=paste(nowotwor, '.clinical', sep=''))

nowotwor.clinical <- nowotwor.clinical[order(nowotwor.clinical$patient.bcr_patient_barcode),]


nowotwor.clinical$patient.days_to_death<- as.numeric(as.character(nowotwor.clinical$patient.days_to_death))
nowotwor.clinical$patient.days_to_last_followup<- as.numeric(as.character(nowotwor.clinical$patient.days_to_last_followup))


nowotwor.geny$status <- ifelse(is.na(nowotwor.clinical$patient.days_to_death),0,1)
nowotwor.geny$time <- ifelse(
  is.na(nowotwor.clinical$patient.days_to_death),
  nowotwor.clinical$patient.days_to_last_followup,
  nowotwor.clinical$patient.days_to_death)

# Sprawdzenie
print(paste('clinical:', nrow(nowotwor.clinical), ', geny:', nrow(nowotwor.geny)))
for(n in 1:nrow(nowotwor.clinical)){
  if( nowotwor.geny$pacjent[n] != nowotwor.clinical$patient.bcr_patient_barcode[n]){
    print(paste('Blad', n))}
 
  
}

if(!dir.exists('Poprawki/Zbiory_v2')){
  dir.create('Poprawki/Zbiory_v2')}
write.table(nowotwor.geny, file=paste('Zbiory/', nowotwor, '.txt', sep=""))

}
