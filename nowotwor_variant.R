nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
nowotwory <- 'BRCA'

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

dir_zbiory <- 'C:/Users/Sebastian/OneDrive/Studia/II stopien/Warsztaty badawcze - pbiecek/'

for(nowotwor in nowotwory){
  
  load(paste(dir_zbiory, 'RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep='')) #Wczytanie pacjentów.
  nowotwor.clinical <- get(paste(nowotwor, '.clinical', sep=""))[ ,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]
  rm(list=paste(nowotwor, '.clinical', sep=''))
  
  load(paste(dir_zbiory, 'RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep='')) #Wczytanie mutacji.
  nowotwor.mutations <- get(paste(nowotwor, '.mutations', sep=""))[ ,c("Hugo_Symbol", "bcr_patient_barcode", "Variant_Classification")]
  rm(list=paste(nowotwor, '.mutations', sep=''))
  
  
  nowotwor.clinical$patient.bcr_patient_barcode <- as.character(nowotwor.clinical$patient.bcr_patient_barcode)
  nowotwor.clinical$status <- ifelse(is.na(nowotwor.clinical$patient.days_to_death),0,1)
  nowotwor.clinical$time <- ifelse(
    is.na(nowotwor.clinical$patient.days_to_death),
    nowotwor.clinical$patient.days_to_last_followup,
    nowotwor.clinical$patient.days_to_death)
  nowotwor.clinical <- nowotwor.clinical[ , c("patient.bcr_patient_barcode", "time", "status")]
  
  nowotwor.mutations$bcr_patient_barcode <- tolower(substr(nowotwor.mutations$bcr_patient_barcode, 1, 12))
  variant_cat <- c("Missense_Mutation" ,"Silent"           , "Frame_Shift_Del" ,  "Frame_Shift_Ins" ,  "In_Frame_Del"    ,  "Nonsense_Mutation",
                   "RNA"          ,     "Splice_Site" ,      "In_Frame_Ins"  ,    "Nonstop_Mutation" )
  variant_id <- 1:10
  variant <- cbind(variant_cat, variant_id)
  colnames(variant) <- c("Kategoria", "Id")
  variant <- data.frame(variant)
  
  tmp <- NULL
  nowotwor.clinical.rep <- sapply(1:3, function(i){
    tmp <- cbind(tmp, rep(nowotwor.clinical[,i], each=10))
    tmp
  })
  nowotwor.clinical.rep <- data.frame(nowotwor.clinical.rep)
  nowotwor.clinical.rep$variant_id <- rep(1:10, 1098)
  
  colnames(nowotwor.clinical.rep) <- c("patient.barcode", "time", "status", "variant_id")
  nowotwor.clinical.rep <- as.data.frame(nowotwor.clinical.rep)
  nowotwor.clinical.rep <- nowotwor.clinical.rep[order(nowotwor.clinical.rep$patient.barcode),]
  
  nowotwor.wynik <- nowotwor.clinical.rep[,c("patient.barcode", "variant_id", "time", "status")]
  nowotwor.wynik <- nowotwor.wynik[order(nowotwor.wynik$patient.barcode, nowotwor.wynik$variant_id),]
  
  k <- 1
  for(gen in geny){
    print(paste(nowotwor, ":", k, "/", length(geny)))
    k <- k + 1
    nowotwor.mutations.gen <- nowotwor.mutations[nowotwor.mutations$Hugo_Symbol == gen, ]
    if(nrow(nowotwor.mutations.gen) == 0){
      nowotwor.wynik_tmp <- cbind(nowotwor.wynik, NA)
    }
    else{
      nowotwor.mutations.gen$variant_id <- 0
      for(v in variant$Kategoria){
        nowotwor.mutations.gen$variant_id[nowotwor.mutations.gen$Variant_Classification == v] <- variant$Id[variant$Kategoria ==v]
      }
      
      nowotwor.merge <- merge(nowotwor.clinical.rep, nowotwor.mutations.gen, by.x = c("patient.barcode", "variant_id"), by.y = c("bcr_patient_barcode", "variant_id"), all.x = T, all.y = F)
      nowotwor.merge <- unique(nowotwor.merge)
      nowotwor.merge <- nowotwor.merge[order(nowotwor.merge$patient.barcode, nowotwor.merge$variant_id),]
      colnames(nowotwor.merge) <- c("patient.barcode", "variant_id", "time", "status", "gen", paste("Variant.", gen, sep=""))
      
      
      nowotwor.wynik_tmp <- cbind(nowotwor.wynik, nowotwor.merge[,ncol(nowotwor.merge)])
    }
    colnames(nowotwor.wynik_tmp) <- c(colnames(nowotwor.wynik), paste("Variant.", gen, sep=""))
    nowotwor.wynik <- nowotwor.wynik_tmp
  }
  
  write.table(nowotwor.wynik, file = paste("Zbiory/", nowotwor, "_variant.txt", sep=""))
  
}

