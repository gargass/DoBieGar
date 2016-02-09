#nowotwory= c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
#             "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
nowotwor= 'BRCA'
dane = read.table(paste('C:/Users/Gosia/Documents/DoBieGar/Zbiory/', nowotwor,'.txt' ,sep=''), h=T)

geny <- read.table('C:/Users/Gosia/Documents/DoBieGar/p_value/lista_interesujacych_genow.txt', h=T)
wspolne = matrix(0, nrow = 535, ncol=535)
geny$x
colnames(wspolne)<-geny$x
rownames(wspolne)<-geny$x


for (i in 1:length(geny))
{
  for (j in i:length(geny))
  {
    if (sum(as.numeric(dane[,i+1]))!=0)
    {
    wspolne[i, j] = sum(as.numeric(dane[dane[, i+1]==1 & dane[, j+1]==1,j+1]))/sum(as.numeric(dane[,i+1]))
    wspolne[j, i]=wspolne[i, j]
    }
    else
    {
      wspolne[i, j]=wspolne[j, i]=0
    }
  }
  
}

write.table(file=paste('Zbiory/wspolne_', nowotwor, '.txt', sep=""), wspolne)
