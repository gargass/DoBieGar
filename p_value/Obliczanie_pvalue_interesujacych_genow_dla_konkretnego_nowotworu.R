library(survival)

setwd('C:/Users/Gosia/Desktop/Biecek')

nowotwor <- 'UCEC'
lista_genow = read.table('C:/Users/Gosia/Desktop/Biecek/Istotne_geny_dla_nowotworu/Interesujace geny/lista_interesujacych_genow.txt', h=T)
#c('BRCA', COADREAD', 'GBM', 'GBMLGG', 'HNSC', 'KIPAN', 'KIRC', 'LGG', 'LUAD', 'LUSC', 
#'OV', 'STES', 'THCA', 'UCEC')

load(paste('RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep='')) #Wczytanie pacjentĂłw.
load(paste('RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep='')) #Wczytanie mutacji.

nowotwor.clinical_sub<-get(paste(nowotwor, '.clinical', sep=""))[,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]
nowotwor.clinical_sub$patient.days_to_death<-as.numeric(as.character(nowotwor.clinical_sub$patient.days_to_death))
nowotwor.clinical_sub$patient.days_to_last_followup<-as.numeric(as.character(nowotwor.clinical_sub$patient.days_to_last_followup))

nowotwor.mutations_sub<-get(paste(nowotwor, '.mutations', sep=""))[!is.na(get(paste(nowotwor, '.mutations', sep=""))$Hugo_Symbol),c("Hugo_Symbol", "Tumor_Sample_Barcode")]
nowotwor.mutations_sub$Tumor_Sample_Barcode<-tolower(substr(nowotwor.mutations_sub$Tumor_Sample_Barcode, 1, 12))
nowotwor.mutations_sub$czy_jest<-1


nowotwor.geny <- lista_genow$x


nowotwor.geny_n<-seq(0,0,length=length(nowotwor.geny))
names(nowotwor.geny_n)<-nowotwor.geny

nowotwor.cli_mut<-nowotwor.clinical_sub
nowotwor.cli_mut$status<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),0,1)
nowotwor.cli_mut$time<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),nowotwor.cli_mut$patient.days_to_last_followup,nowotwor.cli_mut$patient.days_to_death)


nowotwor.P_value<-seq(0,0,length=length(nowotwor.geny))
for(i in 1:length(nowotwor.geny))
  
  {
  
  nowotwor.mutations_sub2<-nowotwor.mutations_sub[nowotwor.mutations_sub$Hugo_Symbol==nowotwor.geny[i],2:3]
  
  if (nrow(nowotwor.mutations_sub2)>0)
{
  colnames(nowotwor.mutations_sub2)<-c('Tumor_Sample_Barcode', nowotwor.geny[i])
  nowotwor.cli_mut_temp<-merge(nowotwor.cli_mut, nowotwor.mutations_sub2, by.y='Tumor_Sample_Barcode', by.x='patient.bcr_patient_barcode', all.x=T, all.y=F)
  nowotwor.cli_mut_temp[,6]<-ifelse(is.na(nowotwor.cli_mut_temp[,6]),0,nowotwor.cli_mut_temp[,6])
  if ((sum(nowotwor.cli_mut_temp[,6])==length(nowotwor.cli_mut_temp[,6])) || (sum(nowotwor.cli_mut_temp[,6])==0))
  {nowotwor.P_value[i]<-1}
  else
  {
    nowotwor.surv<-survdiff(Surv(time, status)~nowotwor.cli_mut_temp[,6], data=nowotwor.cli_mut_temp)
    nowotwor.P_value[i]<-1-pchisq(nowotwor.surv$chisq,1)
  }
}

  else 
    {nowotwor.P_value[i]<-9999}

}
names(nowotwor.P_value)<-nowotwor.geny

 nowotwor.P_value.wypisz<-as.data.frame(round(nowotwor.P_value, 6))
colnames(nowotwor.P_value.wypisz)=c('Pvalue')


if(!dir.exists('P_value_dla_interesujacych_genow')){
  dir.create('P_value_dla_interesujacych_genow')
}


write.table(file=paste('P_value_dla_interesujacych_genow/',  nowotwor, '_pvalue.txt', sep=''), round(nowotwor.P_value.wypisz,6))
