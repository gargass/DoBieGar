library(shiny)
library(survival)
library(survMisc)
library(grid)
library(gridExtra)

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

shinyServer(function(input, output) {
  
  output$opis_wykres <- renderText({
    nowotwor <- input$nowotwory
    gen <- input$geny
    if( nowotwor != 'Wszystkie'){
      p_value <- read.table(paste('p_value/P_value_dla_interesujacych_genow/', 
                                  nowotwor, '_pvalue.txt', sep=""), h=T)
      
      paste("P-value: ", p_value$Pvalue[rownames(p_value) == gen])
    }
  })
  
  output$wykres <- renderPlot({
    nowotwory <- input$nowotwory
    gen <- input$geny
    if(!is.element("Wszystkie", nowotwory)){
      
      n <- length(nowotwory)
     
      p <- lapply(nowotwory, function(nowotwor){
        
        zbior.nowotwor <- read.table(paste('Zbiory/', nowotwor, '.txt', sep=""))
        p_value <- read.table(paste('p_value/P_value_dla_interesujacych_genow/', 
                                    nowotwor, '_pvalue.txt', sep=""), h=T)
        nowotwor_gen.fit <- survfit(Surv(time, status) ~ zbior.nowotwor[,gen], 
                                    data=zbior.nowotwor)
        survMisc::autoplot(nowotwor_gen.fit,
                           xLab = paste('Time, P-value:', p_value$Pvalue[rownames(p_value) == gen]),
                           legLabs = c("status = 0","status = 1"),
                           title=paste('Krzywa przeżycia dla genu ', gen, '\n w nowotworze ', nowotwor, sep=""))$plot
      })

      if(n <= 4){
        ncol <- 2
        nrow <- 2
      }
      else{
        ncol <- ceiling(sqrt(n))
        nrow <- ceiling(n/ceiling(sqrt(n)))
      }
      
      marrangeGrob(p, ncol = ncol, nrow = nrow)
      
    }}, height = 1000, width = 1000)
  
  output$opis_geny <- renderText({
    "Tabela przedstawiająca 10 najistotniejszych genów, na których 
    wystąpiła mutacja w nowotworze. \n"
  })
  
  output$geny <- renderTable({
    nowotwor <- input$nowotwory
    if (input$nowotwory!='Wszystkie')
    {
      p_value_tabela <-read.table('p_value/p_value_NA.txt', h=T)
      p <- matrix(1, nrow=10, ncol=1)
      for (nowotworr in nowotwor)
      {
        p2 <- p_value_tabela[order(p_value_tabela[,nowotworr]), ][1:10, c("gen", nowotworr)]
        colnames(p2) <- c(paste("Marker ", nowotworr), paste("P-value ", nowotworr))
        p <- cbind(p, p2)
      }
      p <- p[, -c(1)]
      rownames(p) <- NULL
      print(p)
    }}, digits = 5)
  
  output$opis_geny_wspol <- renderText({
    "Geny współistniejące z wybranymi:"
  })
  
  output$geny_wspolne<-renderTable({
    p <- matrix(1, nrow=3, ncol=1)
    gen <- input$geny
    nowotwory <- input$nowotwory
    for (nowotwor in nowotwory){
      dane <- read.table(paste('Zbiory/', nowotwor, '.txt', sep=''), h=T)
      z <- numeric((ncol(dane)-2))
      for (i in 3:(ncol(dane)-2)){
        z[i] <- sum(dane[which(dane[, gen]==1), i])
        }
      a <- sum(dane[which(dane[, gen]==1), gen])
      x <- z[order(z)][(length(z)-1):(length(z)-3)]
      x2 <- c(which(z==x[1]), which(z==x[2]),which(z==x[3]))
      x2 <- x2[1:3]
      
      x3 <- as.matrix(colnames(dane)[c(x2)])
      x4 <- as.matrix(round(z [ c(x2)]/a, 2))
      colnames(x3) <- c(paste("Marker ", nowotwor))
      colnames(x4) <- c(paste("Correlation ", nowotwor))
      p <- cbind(p, x3, x4)
      }
    p <- p[, -c(1)]
    rownames(p) <- NULL
    print(p)
    })
  })