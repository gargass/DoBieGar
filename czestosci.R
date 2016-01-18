nowotwory= c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
"UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
nowotwor= 'BRCA'
dane = read.table(paste('Zbiory/', nowotwor,'.txt' ,sep=''), h=T)


czestosci = matrix(nrow = (ncol(dane)-3), ncol=15)
czestosci[,1]=colnames(dane)[2:(ncol(dane)-2)]

for (j in 1:length(nowotwory))
{
  
  dane = read.table(paste('Zbiory/', nowotwory[j],'.txt' ,sep=''), h=T)
  n = nrow(dane)
for (i in 1:297)
{
  czestosci[i, j+1]=round(sum(dane[, i+1])/n, 3)
  
}


}

colnames(czestosci)=c('gen', nowotwory)


write.table(file=paste('czestosci.txt', sep=""), czestosci)

