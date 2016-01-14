library(shiny)
library(survival)
library(survMisc)

shinyServer(function(input, output) {
  
  
  setwd<-'C:/Users/Gosia/Desktop/Biecek'
  setwd(setwd)
  dane <- read.table(paste(setwd,'/STES','/p_value.txt', sep=''))
  
  geny <- as.factor(rownames(dane))
  
  
  # tabela
  output$tabela <- renderTable({
    #     serialeIMDB[serialeIMDB$serial == input$serial, c("serial", "nazwa", "sezon", "ocena", "glosow")]
  })
  
  # wykres
  output$wykres <- renderPlot({
    nowotwor<-input$nowotwory
    if (nowotwor != "Wszystkie"){
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
      nowotwor.fit<-survfit(Surv(time, status)~nowotwor.cli_mut_temp[,6], data=nowotwor.cli_mut_temp)
      
      #     tmp <- serialeIMDB[serialeIMDB$serial == input$serial,]
      #     ggplot(tmp, aes(x=id, y=ocena)) + geom_point() + geom_smooth(se=FALSE)
      #plot(nowotwor.fit, 
      #xlab="Dni", ylab="P-stwo przeżycia")
      #legend("bottomleft", c(paste(names(nowotwor.cli_mut_temp[6]), '=0'),
      #paste(names(nowotwor.cli_mut_temp[6]), '=1')), 
      #col=c(1,2), lty=1)
      p = autoplot(nowotwor.fit, title=paste('Krzywa przeżycia dla genu ', input$geny, '\n w nowotworze ', input$nowotwory, sep=""), legLabs=c("status=0", "status=1"))$plot
      
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
