
setwd('C:/Users/Gosia/Desktop/Biecek')
nowotwor <- 'LUSC'

geny = read.table("C:/Users/Gosia/Desktop/Biecek - kody/Pliki z genami istotnymi dla kazdego nowotworu + lista genow, ktore bierzemy pod uwage/lista_interesujacych_genow.txt", h=T)
gen = geny$x[1]

dane = read.table(paste('Zbiory/', nowotwor, '/', nowotwor, '_', gen, '_cli_mut.txt', sep=""))
dane = dane[order(dane$patient.bcr_patient_barcode),]
dane  = dane[!duplicated(dane), ]
for (i in 2:length(geny$x))
{
  gen = geny$x[i]
  dane1 = read.table(paste('Zbiory/', nowotwor, '/', nowotwor, '_', gen, '_cli_mut.txt', sep=""))
  dane1 = dane1[order(dane1$patient.bcr_patient_barcode),]
  dane1  = dane1[!duplicated(dane1), ]
  dane = cbind(dane, dane1[, 2])
  
}
colnames(dane)=c('pacjent', as.character(geny$x))

write.table(dane, file = paste('Zbiory/', nowotwor, '.txt', sep=""))
