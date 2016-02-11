nowotwory <- c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
               "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")


for(nowotwor in nowotwory){
  dane <- read.table(paste('Zbiory/', nowotwor,'.txt' ,sep=''), h=T)
  
  geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
  wspolne <- data.frame(matrix(0, nrow = nrow(geny), ncol=nrow(geny)))
  
  colnames(wspolne)<-geny$x
  rownames(wspolne)<-geny$x
  
  
  for (i in 1:nrow(geny))
  {
    gen_i <- as.character(geny$x[i])
    for (j in i:nrow(geny))
    {
      gen_j <- as.character(geny$x[j])
      
      n.gen_j <- sum(as.numeric(dane[,gen_j] == 1 & dane[, gen_i] ==1))
        wspolne[i, j] <- n.gen_j
        wspolne[j, i] <- n.gen_j 
      
      
    }
    print(paste(nowotwor, ':',i, '/', nrow(geny)))
    
    
  }
  
  write.table(file=paste('Zbiory/wspolne_', nowotwor, '_licznosci.txt', sep=""), wspolne)
}
