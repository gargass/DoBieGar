library(shiny)
library(survival)
library(survMisc)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)

nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

for(nowotwor in nowotwory){
  assign(paste('zbior.', nowotwor, sep=""), read.table(paste('Zbiory/', nowotwor, '.txt', sep="")))
}

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

p_value_tabela <-read.table('p_value/p_value_NA.txt', h=T)
czestosci<-read.table('czestosci.txt', h=T)

for(nowotwor in nowotwory){
  assign(paste('p_value.', nowotwor, sep=""), read.table(paste('p_value/P_value_dla_interesujacych_genow/', 
                                                               nowotwor, '_pvalue.txt', sep=""), h=T))
}

for(nowotwor in nowotwory){
  assign(paste(nowotwor, '_variant', sep=""), read.table(paste('Zbiory/', nowotwor, '_variant.txt', sep="")))
}

najczestsze <- read.table("najistotniejsze_geny.txt", h=T)

shinyServer(function(input, output) {
  
  output$opis_krzywe <- renderText({
    nowotwor <- input$nowotwory
    gen <- input$geny

    print('In the figures below we can see Kaplan-Meier curves for 
          a given gene and the given tumors. The survival curves 
          are estimated for the two groups of patients: 
          the first one refers to the patients with a mutation of a given gene 
          and the second one is the group of patients without any mutation 
          of this gene.')
  })
  

  output$wykres <- renderPlot({
    nowotwory <- input$nowotwory
    gen <- input$geny
  
      
      n <- length(nowotwory)
     
      p <- lapply(nowotwory, function(nowotwor){
        
        nowotwor_gen.fit <- survfit(Surv(time, status) ~ get(paste('zbior.', nowotwor, sep=""))[,gen], 
                                    data=get(paste('zbior.', nowotwor, sep="")))
        survMisc::autoplot(nowotwor_gen.fit,
                           xLab = paste('Time, P-value:', 
                                        get(paste('p_value.', nowotwor, sep=""))$Pvalue[rownames(get(paste('p_value.', nowotwor, sep=""))) == gen]), 
                           legLabs = c("status = 0","status = 1"),
                           title=paste('Kaplan–Meier estimator for ', gen, '\n in ', nowotwor, " cancer", sep=""))$plot
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
      
    }, height = 1000, width = 1000)
  
  output$opis_geny <- renderText({
    "The table shows the 10 most significant genes in which 
    a mutation occurred in cancer. \n"
  })
  
  output$geny <- renderTable({
    nowotwor <- input$nowotwory
   
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
    }, digits = 5)
  
  output$opis_geny_wspol <- renderText({
    "The following table depicts three genes whose mutations appear most frequently with the mutation of a given gene in the given tumors."
  })
  
  output$geny_wspolne<-renderTable({
    p <- matrix(1, nrow=3, ncol=1)
    gen <- input$geny
    nowotwory <- input$nowotwory
    for (nowotwor in nowotwory){
      dane <- get(paste('zbior.', nowotwor, sep=""))
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
  
  output$heatmap_pvalue <- renderPlot({
    melted_dane <- melt(p_value_tabela[which(p_value_tabela$gen %in% najczestsze$x[1:50]), ])
    
    base_size <- 12
    
    ggplot(data = melted_dane, aes(x=variable, y=gen, fill=value)) + 
      geom_tile() + theme_grey(base_size = base_size) + labs(x = "",y = "") + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(), 
                                                 axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))+
      geom_tile(aes(fill = value), colour = "white")
    
  })
  
  output$heatmap_czestosc <- renderPlot({
    melted_dane <- melt(czestosci[which(czestosci$gen %in% najczestsze$x[1:50]), ])
    
    base_size <- 12
    
    
    ggplot(data = melted_dane, aes(x=variable, y=gen, fill=value)) + 
      geom_tile() + theme_grey(base_size = base_size) + labs(x = "",y = "") + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(), 
                                                 axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))+
      geom_tile(aes(fill = value), colour = "white") + 
      scale_fill_gradient(low = "white", high = "steelblue")
    
  })
  
  output$boxplot_variant <- renderPlot({
      nowotwory <- input$nowotwory
      gen <- input$geny
      
      n <- length(nowotwory)
      
      nowotwory.gen.missense <- NULL
      nowotwory.gen.nonsense <- NULL
      
      for(nowotwor in nowotwory){
        nowotwor_variant <- get(paste(nowotwor, '_variant', sep=""))
        nowotwor_variant <- nowotwor_variant[, c("patient.barcode", "time", "status", paste('Variant.', gen, sep=""))]
        nowotwor_variant$nowotwor <- nowotwor
        colnames(nowotwor_variant) <- c("patient.barcode", "time", "status", "Variant", "nowotwor")
        nowotwory.gen.missense <- rbind(nowotwory.gen.missense, nowotwor_variant[!is.na(nowotwor_variant$Variant) & nowotwor_variant$Variant== "Missense_Mutation", ])
        nowotwory.gen.nonsense <- rbind(nowotwory.gen.nonsense, nowotwor_variant[!is.na(nowotwor_variant$Variant) & nowotwor_variant$Variant == "Nonsense_Mutation", ])
        
      }
        
      #ggplot(nowotwory.gen.missense, aes(x=nowotwor, y=time)) + geom_boxplot()
      par(mfrow = c(1,2))
      boxplot(time ~ nowotwor, data = nowotwory.gen.missense, main = "Missense Mutation")
      boxplot(time ~ nowotwor, data = nowotwory.gen.nonsense, main = "Nonsense Mutation")
      
      
     
      
  }, height = 400, width = 800)
  })