library(shiny)
library(survival)
library(survMisc)

#zal zal zal
#TRAAAAAAAALALLALLAa

#nowotwory = c('COADREAD', 'GBM', 'GBMLGG', 'HNSC', 'KIPAN', 'KIRC', 'LGG', 'LUAD', 'LUSC',
#               'OV', 'STES', 'THCA', 'UCEC')
nowotwory <- list("GBMLGG", "BRCA", "KIPAN")

dane = NULL
for (nowotwor in nowotwory)
{
  dane1 = read.table(paste('pvalue/', nowotwor, '_pvalue.txt', sep=''), h=T)
  dane1$gen = rownames(dane1)
  dane1$nowotwor = nowotwor
  dane = rbind(dane, dane1)
}
geny<-unique(dane$gen)


shinyServer(function(input, output) {
  
  
  #   setwd<-'C:/Users/Gosia/Desktop/Biecek'
  #   setwd(setwd)
  #  dane <- read.table(paste(setwd,'/STES','/p_value.txt', sep=''))
  
  #   dane <- read.table('KIPAN/KIPAN_p_value.txt', sep='')
  #   geny <- as.factor(rownames(dane))
  # tabela
  output$tabela <- renderTable({
    gen <- input$geny
    
    #     serialeIMDB[serialeIMDB$serial == input$serial, c("serial", "nazwa", "sezon", "ocena", "glosow")]
    if (input$nowotwory =='Wszystkie'){
      
      poziomy_genu <- NULL
      for(n in 1:length(nowotwory)){
        zestaw.nowotwor_gen <- read.table(paste('Zbiory/', nowotwory[n], '/', nowotwory[n], '_', gen, '_cli_mut.txt', sep=""))
        poziomy_genu <- c(poziomy_genu, as.character(unique(zestaw.nowotwor_gen$Variant_Classification)))
      }
      poziomy_genu <- unique(poziomy_genu)
      
      tabela <- matrix(0, nrow = length(nowotwory), ncol=3 + length(poziomy_genu))
      rownames(tabela) <- nowotwory
      colnames(tabela) <- c('Liczność genu w zbiorze', 'Liczność zbioru', 'Frakcja wystąpień', poziomy_genu)
      
      for(n in 1:length(nowotwory)){
        zestaw.nowotwor_gen <- read.table(paste('Zbiory/', nowotwory[n], '/', nowotwory[n], '_', gen, '_cli_mut.txt', sep=""))
        
        tabela[n,1] <- sum(zestaw.nowotwor_gen[,ncol(zestaw.nowotwor_gen)])
        tabela[n,2] <- nrow(zestaw.nowotwor_gen)
        tabela[n,3] <- tabela[n,1] / tabela[n,2]
        for(poz in 1:length(poziomy_genu)){
          
          #tabela[n, 3 + poz] <- nrow(zestaw.nowotwor_gen[as.character(zestaw.nowotwor_gen$Variant_Classification) == poziomy_genu[poz] & zestaw.nowotwor_gen$status==1,])/tabela[n,1]
          tabela[n, 3 + poz] <- nrow(zestaw.nowotwor_gen[as.character(zestaw.nowotwor_gen$Variant_Classification) == poziomy_genu[poz] & zestaw.nowotwor_gen[,ncol(zestaw.nowotwor_gen)]==1,])/tabela[n,1]
        }
        
        
      }
      print(tabela)
      
    }
    
    
  })
  
  output$tabela_kilka <- renderTable({
    gen <- input$geny
    tabela <- matrix(0, nrow=length(nowotwory), ncol = length(nowotwory))
    rownames(tabela) <- nowotwory
    colnames(tabela) <- nowotwory
    
    for(t1 in 1:length(nowotwory)){
      p_value_t1 <- read.table(paste('pvalue/', nowotwory[t1], '_pvalue.txt', sep=''), h=T)
      for(t2 in 1:length(nowotwory)){
        p_value_t2 <- read.table(paste('pvalue/', nowotwory[t2], '_pvalue.txt', sep=''), h=T)
        tabela[t1,t2] <- (is.element(gen, rownames(p_value_t1)) & is.element(gen, rownames(p_value_t2)))
      }
    }
    
    print(tabela)
    
  })
  
  # wykres
  output$wykres <- renderPlot({
    nowotwor<-input$nowotwory
    gen <- input$geny
    if (nowotwor != "Wszystkie"){
      
      zestaw.nowotwor_gen <- read.table(paste('Zbiory/', nowotwor, '/', nowotwor, '_', gen, '_cli_mut.txt', sep=""))
      nowotwor_gen.fit<-survfit(Surv(time, status) ~ zestaw.nowotwor_gen[,7], data=zestaw.nowotwor_gen)
      
      #     tmp <- serialeIMDB[serialeIMDB$serial == input$serial,]
      #     ggplot(tmp, aes(x=id, y=ocena)) + geom_point() + geom_smooth(se=FALSE)
      #plot(nowotwor.fit,
      #xlab="Dni", ylab="P-stwo przeżycia")
      #legend("bottomleft", c(paste(names(nowotwor.cli_mut_temp[6]), '=0'),
      #paste(names(nowotwor.cli_mut_temp[6]), '=1')),
      #col=c(1,2), lty=1)
      p = autoplot(nowotwor_gen.fit, title=paste('Krzywa przeżycia dla genu ', input$geny, '\n w nowotworze ', input$nowotwory, sep=""), legLabs=c("status=0", "status=1"))$plot
      
      print(p)
    }
    
  })
  
  # podsumowanie
  output$podsumowanie <- renderPrint({
    nowotwor=input$nowotwory
    load(paste('RTCGA.clinical-master/data/', nowotwor, '.clinical.rda', sep=''))
    load(paste('RTCGA.mutations-master/data/', nowotwor, '.mutations.rda', sep=''))
    nowotwor.clinical<-get(paste(nowotwor, '.clinical', sep=''))
    nowotwor.clinical_sub<-nowotwor.clinical[,c("patient.bcr_patient_barcode", "patient.days_to_death", "patient.days_to_last_followup")]
    
    nowotwor.mutations<-get(paste(nowotwor, '.mutations', sep=''))
    nowotwor.mutations_sub<-nowotwor.mutations[!is.na(nowotwor.mutations$Hugo_Symbol),c("Hugo_Symbol", "Tumor_Sample_Barcode")]
    nowotwor.mutations_sub$Tumor_Sample_Barcode<-tolower(substr(nowotwor.mutations_sub$Tumor_Sample_Barcode, 1, 12))
    nowotwor.mutations_sub$czy_jest<-1
    
    nowotwor.cli_mut<-nowotwor.clinical_sub
    nowotwor.cli_mut$status<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),0,1)
    nowotwor.cli_mut$time<-ifelse(is.na(nowotwor.cli_mut$patient.days_to_death),nowotwor.cli_mut$patient.days_to_last_followup,nowotwor.cli_mut$patient.days_to_death)
    
    nowotwor.mutations_sub2<-nowotwor.mutations_sub[nowotwor.mutations_sub$Hugo_Symbol==input$geny,2:3]
    colnames(nowotwor.mutations_sub2)<-c('Tumor_Sample_Barcode', input$geny)
    nowotwor.cli_mut_temp<-merge(nowotwor.cli_mut, nowotwor.mutations_sub2, by.y='Tumor_Sample_Barcode', by.x='patient.bcr_patient_barcode', all.x=T, all.y=F)
    nowotwor.cli_mut_temp[,6]<-ifelse(is.na(nowotwor.cli_mut_temp[,6]),0,nowotwor.cli_mut_temp[,6])
    
    
    x=nrow(nowotwor.clinical_sub)
    y=sum(nowotwor.cli_mut_temp[,6])
    z=y/x
    w=cbind(x, y, z)
    colnames(w)=c("Liczba pacjentów", "Liczba pacjentów z mutacją genu", "Procent")
    print(w)                        
    #tmp <- serialeIMDB[serialeIMDB$serial == input$serial, c("ocena", "glosow")]
    #summary(tmp)
  })
  
})
