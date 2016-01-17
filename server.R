library(shiny)
library(survival)
library(survMisc)
library(grid)
library(gridExtra)

#zal zal zal zal zal
#TRAAAAAAAALALLALLAa

#nowotwory = c('COADREAD', 'GBM', 'GBMLGG', 'HNSC', 'KIPAN', 'KIRC', 'LGG', 'LUAD', 'LUSC',
#               'OV', 'STES', 'THCA', 'UCEC')
# nowotwory <- list("GBMLGG", "BRCA", "KIPAN")
# 
# dane = NULL
# for (nowotwor in nowotwory)
# {
#   dane1 = read.table(paste('pvalue/', nowotwor, '_pvalue.txt', sep=''), h=T)
#   dane1$gen = rownames(dane1)
#   dane1$nowotwor = nowotwor
#   dane = rbind(dane, dane1)
# }
# geny<-unique(dane$gen)
geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

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
  output$text1 <- renderText({ 
    nowotwor<-input$nowotwory
    gen <- input$geny
    if( nowotwor != 'Wszystkie'){
      p_value <-
        read.table(
          paste('p_value/P_value_dla_interesujacych_genow/', nowotwor, '_pvalue.txt', sep=""), 
          h=T)
      
      paste("P-value: ", p_value$Pvalue[rownames(p_value) == gen])
    }
  })
  
  output$wykres <- renderPlot({
    nowotwory <- input$nowotwory
    gen <- input$geny
    #if (nowotwor != "Wszystkie"){
    if(!is.element("Wszystkie", nowotwory)){
      n <- length(nowotwory)
      #print(n)
      #par(mfrow=c(floor(sqrt(n)), ceiling(sqrt(n))))
#       grobs <- gList()
#       grobs_i <- 1
#       for(nowotwor in nowotwory){
#         print(nowotwor)
#         zbior.nowotwor<- read.table(paste('Zbiory/', nowotwor, '.txt', sep=""))
#         nowotwor_gen.fit<-survfit(Surv(time, status) ~ zbior.nowotwor[,gen], data=zbior.nowotwor)
#         
#         
#         grobs[grobs_i] <- grob(autoplot(
#           nowotwor_gen.fit, 
#           title=paste('Krzywa przeżycia dla genu ', gen, '\n w nowotworze ', nowotwor, sep=""), 
#           legLabs=c("status=0", "status=1"))$plot)
#         grobs_i = grobs_i + 1
#       }
#       
#       arrangeGrob(grobs)
        p <- lapply(nowotwory, function(nowotwor) {
          zbior.nowotwor<- read.table(paste('Zbiory/', nowotwor, '.txt', sep=""))
          nowotwor_gen.fit<-survfit(Surv(time, status) ~ zbior.nowotwor[,gen], data=zbior.nowotwor)
          
          survMisc::autoplot(nowotwor_gen.fit, legLabs = c("lower","higher"))$plot 
          })
        
        marrangeGrob(p, nrow = floor(sqrt(n)), ncol=ceiling(sqrt(n)))
      
    }
      
    
    
    
    
  }, height = 1000, width = 1000)
  
  # podsumowanie
  output$geny <- renderTable({
    nowotwor = input$nowotwory
    if (input$nowotwory!='Wszystkie')
    {
      
      # nowotwor=c('BRCA', 'LGG')
      p_value_tabela <-read.table('p_value/p_value_NA.txt', h=T)
      p = matrix(1, nrow=20, ncol=1)
      for (nowotworr in nowotwor)
      {
      
        p2=p_value_tabela[order(p_value_tabela[,nowotworr]), ][1:20, c("gen", nowotworr)]
        p = cbind(p, p2)
        
      }
      
      #p1 =as.table(p1)
      p = p[, -c(1)]
      p <- as.data.frame(p)
      #rownames(p)=NULL
      
      print(p)
      #p=p[, -c(1, 2)]

    }
    
  }, digits = 5)
  
})
