nowotwory= c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
             "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
varianty <- c("Missense_Mutation", "Silent", "Frame_Shift_Del", "Frame_Shift_Ins", 
              "In_Frame_Del", "Nonsense_Mutation", "RNA", "Splice_Site",  
              "In_Frame_Ins", "Nonstop_Mutation")

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

for(nowotwor in nowotwory){
  assign(paste(nowotwor, '_variant', sep=""), read.table(paste('Zbiory/', nowotwor, '_variant.txt', sep="")))
}


czestosci = matrix(nrow = length(varianty)*length(nowotwory), ncol=length(geny)+2)

for (i in 1:length(nowotwory)){
  dane <- get(paste(nowotwory[i], '_variant', sep=""))
  n = nrow(dane)/10
  for (j in 1:length(varianty)){
    czestosci[(i - 1)*length(varianty) + j, 1] <- nowotwory[i]
    czestosci[(i - 1)*length(varianty) + j, 2] <- varianty[j]
    for (k in 1:length(geny)){
      czestosci[(i - 1)*length(varianty) + j, k+2] <- 
        round(nrow(dane[!is.na(dane[,k+4]) & dane[,k+4]==varianty[j], ])/n, 3)
    }
  }
}

colnames(czestosci)=c('nowotwor', 'variant', geny)

write.table(file=paste('czestosci_variant.txt', sep=""), czestosci)
