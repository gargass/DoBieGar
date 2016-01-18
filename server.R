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
  
  
  

  
  # wykres
  output$opis_wykres <- renderText({ 
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
     
      p <- lapply(nowotwory, function(nowotwor) {
        zbior.nowotwor<- read.table(paste('Zbiory/', nowotwor, '.txt', sep=""))
        p_value <-
          read.table(
            paste('p_value/P_value_dla_interesujacych_genow/', nowotwor, '_pvalue.txt', sep=""), 
            h=T)
        nowotwor_gen.fit<-survfit(Surv(time, status) ~ zbior.nowotwor[,gen], data=zbior.nowotwor)
        survMisc::autoplot(nowotwor_gen.fit,
                           xLab = paste('Time, P-value:', p_value$Pvalue[rownames(p_value) == gen]),
                           legLabs = c("status = 0","status = 1"),
                           title=paste('Krzywa przeżycia dla genu ', gen, '\n w nowotworze ', nowotwor, sep=""))$plot

        

       
      })

      if(n<=4){
        ncol = 2
        nrow = 2
      }
      else{
        ncol = ceiling(sqrt(n))
        nrow = ceiling(n/ceiling(sqrt(n)))
      }
      
      marrangeGrob(p, ncol = ncol, nrow = nrow)
      
    }
    
    
    
    
    
  }, height = 1000, width = 1000)
  
  # podsumowanie
  output$opis_geny <- renderText({
    "Tabela przedstawiająca 10 najistotniejszych genów, na których wystąpiła mutacja w nowotworze \n"
  })
  output$geny <- renderTable({
    nowotwor = input$nowotwory
    if (input$nowotwory!='Wszystkie')
    {
      
      # nowotwor=c('BRCA', 'LGG')
      p_value_tabela <-read.table('p_value/p_value_NA.txt', h=T)
      p = matrix(1, nrow=10, ncol=1)
      for (nowotworr in nowotwor)
      {
        
        p2=p_value_tabela[order(p_value_tabela[,nowotworr]), ][1:10, c("gen", nowotworr)]
        colnames(p2) <- c(paste("Marker ", nowotworr), paste("P-value ", nowotworr))
        p = cbind(p, p2)
        
      }
      
      #p1 =as.table(p1)
      p = p[, -c(1)]
      rownames(p)=NULL
      print(p)
      #p=p[, -c(1, 2)]
      
    }
    
  }, digits = 5)
  
  output$opis_geny_wspol <- renderText({
    "Geny współistniejące z wybranymi"
  })
  output$geny_wspolne<-renderTable({
    p = matrix(1, nrow=3, ncol=1)
 gen = input$geny
  nowotwory=input$nowotwory
 for (nowotwor in nowotwory)
 {
  dane = read.table(paste('Zbiory/', nowotwor, '.txt', sep=''), h=T)
  
  z = numeric((ncol(dane)-2))
  for (i in 3:(ncol(dane)-2))
  {
    
    z[i]=sum(dane[which(dane[, gen]==1), i])
    
  }
  a = sum(dane[which(dane[, gen]==1), gen])
  x = z[order(z)][(length(z)-1):(length(z)-3)]
  x2 = c(which(z==x[1]), which(z==x[2]),which(z==x[3]))
  x2 = x2[1:3]
  
  x3 = as.matrix(colnames(dane)[c(x2)])
  x4 = as.matrix(round(z [ c(x2)]/a, 2))
  colnames(x3) <- c(paste("Marker ", nowotwor))
  colnames(x4) <- c(paste("Correlation ", nowotwor))
  p =cbind(p, x3, x4)
 }
 p = p[, -c(1)]
 rownames(p)=NULL
 

 print(p)
  
  })
  
#   output$tabela <- renderTable({
#     gen <- input$geny
#     
#     #     serialeIMDB[serialeIMDB$serial == input$serial, c("serial", "nazwa", "sezon", "ocena", "glosow")]
#     if (input$nowotwory =='Wszystkie'){
#       
#       poziomy_genu <- NULL
#       for(n in 1:length(nowotwory)){
#         zestaw.nowotwor_gen <- read.table(paste('Zbiory/', nowotwory[n], '/', nowotwory[n], '_', gen, '_cli_mut.txt', sep=""))
#         poziomy_genu <- c(poziomy_genu, as.character(unique(zestaw.nowotwor_gen$Variant_Classification)))
#       }
#       poziomy_genu <- unique(poziomy_genu)
#       
#       tabela <- matrix(0, nrow = length(nowotwory), ncol=3 + length(poziomy_genu))
#       rownames(tabela) <- nowotwory
#       colnames(tabela) <- c('Liczność genu w zbiorze', 'Liczność zbioru', 'Frakcja wystąpień', poziomy_genu)
#       
#       for(n in 1:length(nowotwory)){
#         zestaw.nowotwor_gen <- read.table(paste('Zbiory/', nowotwory[n], '/', nowotwory[n], '_', gen, '_cli_mut.txt', sep=""))
#         
#         tabela[n,1] <- sum(zestaw.nowotwor_gen[,ncol(zestaw.nowotwor_gen)])
#         tabela[n,2] <- nrow(zestaw.nowotwor_gen)
#         tabela[n,3] <- tabela[n,1] / tabela[n,2]
#         for(poz in 1:length(poziomy_genu)){
#           
#           #tabela[n, 3 + poz] <- nrow(zestaw.nowotwor_gen[as.character(zestaw.nowotwor_gen$Variant_Classification) == poziomy_genu[poz] & zestaw.nowotwor_gen$status==1,])/tabela[n,1]
#           
#           tabela[n, 3 + poz] <- nrow(zestaw.nowotwor_gen[as.character(zestaw.nowotwor_gen$Variant_Classification) == poziomy_genu[poz] & zestaw.nowotwor_gen[,ncol(zestaw.nowotwor_gen)]==1,])/tabela[n,1]
#         }
#         
#         
#       }
#       print(tabela)
#       
#     }
#     
#     
#   })
#   
#   output$tabela_kilka <- renderTable({
#     gen <- input$geny
#     tabela <- matrix(0, nrow=length(nowotwory), ncol = length(nowotwory))
#     rownames(tabela) <- nowotwory
#     colnames(tabela) <- nowotwory
#     
#     for(t1 in 1:length(nowotwory)){
#       p_value_t1 <- read.table(paste('pvalue/', nowotwory[t1], '_pvalue.txt', sep=''), h=T)
#       for(t2 in 1:length(nowotwory)){
#         p_value_t2 <- read.table(paste('pvalue/', nowotwory[t2], '_pvalue.txt', sep=''), h=T)
#         tabela[t1,t2] <- (is.element(gen, rownames(p_value_t1)) & is.element(gen, rownames(p_value_t2)))
#       }
#     }
#     
#     print(tabela)
#     
#   }) 
  
  
})