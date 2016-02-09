nowotwory= c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
             "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
nowotwor= 'BRCA'
dane = read.table(paste('Zbiory/', nowotwor,'.txt' ,sep=''), h=T)


licznosci = matrix(nrow = (ncol(dane)-3), ncol=15)
licznosci[,1]=colnames(dane)[2:(ncol(dane)-2)]

for (j in 1:length(nowotwory))
{
  
  dane = read.table(paste('Zbiory/', nowotwory[j],'.txt' ,sep=''), h=T)
  n = nrow(dane)
  for (i in 1:535)
  {
    licznosci[i, j+1]=sum(dane[, i+1])
    
  }
  
  
}

colnames(licznosci)=c('gen', nowotwory)


write.table(file=paste('licznosci.txt', sep=""), licznosci)
