nowotwor <- 'UCEC'


geny = read.table("p_value/lista_interesujacych_genow.txt", h=T)


load(paste('RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep='')) #Wczytanie pacjentów.
load(paste('RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep='')) #Wczytanie mutacji.

nowotwor.clinical_sub<-get(paste(nowotwor, '.clinical', sep=""))[,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]

nowotwor.mutations_sub<-get(paste(nowotwor, '.mutations', sep=""))[!is.na(get(paste(nowotwor, '.mutations', sep=""))$Hugo_Symbol),c("Hugo_Symbol", "Tumor_Sample_Barcode")]
nowotwor.mutations_sub$Tumor_Sample_Barcode<-tolo/home/marley/Documents/Warsztaty - Biecek/wer(substr(nowotwor.mutations_sub$Tumor_Sample_Barcode, 1, 12))

nowotwor.mutations_sub$czy_jest<-1

nowotwor.cli_mut<-nowotwor.clinical_sub
nowotwor.cli_mut$status<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),0,1)
nowotwor.cli_mut$time<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),nowotwor.cli_mut$patient.days_to_last_followup,nowotwor.cli_mut$patient.days_to_death)

for (i in 1:length(geny$x))
{

  gen = geny$x[i]
nowotwor.mutations_sub2<-nowotwor.mutations_sub[nowotwor.mutations_sub$Hugo_Symbol==gen,2:3]
colnames(nowotwor.mutations_sub2) <- c('Tumor_Sample_Barcode', as.character(gen))
nowotwor.cli_mut_temp <- merge(nowotwor.cli_mut, nowotwor.mutations_sub2, by.y='Tumor_Sample_Barcode', by.x='patient.bcr_patient_barcode', all.x=T, all.y=F)
k <- ncol(nowotwor.cli_mut_temp)
nowotwor.cli_mut_temp[,k]<-ifelse(is.na(nowotwor.cli_mut_temp[,k]),0,nowotwor.cli_mut_temp[,k])


#gdy inny niż pierwszy nowotwor
nowotwor.cli_mut_temp= nowotwor.cli_mut_temp[, c(1, 6)]

if(!dir.exists('Zbiory')){
  dir.create('Zbiory')}

if(!dir.exists(paste('Zbiory/', nowotwor, sep=""))){
  dir.create(paste('Zbiory/', nowotwor, sep=""))}

write.table(nowotwor.cli_mut_temp, file = paste('Zbiory/', nowotwor, '/', nowotwor, '_', gen, '_cli_mut.txt', sep=""))
 }
#}

#}