library(shiny)
library(survival)
# library(survMisc)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
#install.packages('ztable')
#install.packages('pheatmap')
library(pheatmap)
library(DT)
##### Załadowanie zbiorów #####
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

czestosci_variant <- read.table("czestosci_variant.txt", h=T)

#####
shinyServer(function(input, output) {

  
  output$opis_krzywe <- renderText({
    print('In the figures below we can see Kaplan-Meier curves for 
          a given gene and the given tumors. The survival curves 
          are estimated for the two groups of patients: 
          the first one refers to the patients with a mutation of a given gene 
          and the second one is the group of patients without any mutation 
          of this gene.')
  })
  

  output$survcurves_yesno <- renderPlot({
    validate(
      need(input$nowotwory != "", "Please select a cancer!")
    )
    nowotwory <- input$nowotwory
    gen <- input$geny
      
      max_time <- 0
      for(nowotwor in nowotwory){
        zbior <- get(paste('zbior.', nowotwor, sep=""))
        time <- as.numeric(as.character(zbior$time))
        time <- max(time, na.rm = TRUE)
        max_time <- ifelse(time>max_time, time, max_time)
      }    
      
      n <- length(nowotwory)
     
      p <- lapply(nowotwory, function(nowotwor){
        pvalue = p_value_tabela[p_value_tabela$gen == gen, nowotwor]
        pvalue <- signif(pvalue, 3)
        
        nowotwor_gen.fit <- survfit(Surv(as.numeric(as.character(time)), status) ~ get(paste('zbior.', nowotwor, sep=""))[,gen], 
                                    data=get(paste('zbior.', nowotwor, sep="")))
        
        if (sum(get(paste('zbior.', nowotwor, sep=""))[,gen])==0){
          mutation <- "No Mutation"
        }
        else{
          mutation <- c("No Mutation", "Mutation")
        }
      
        
        survMisc::autoplot(nowotwor_gen.fit,
                           legLabs = mutation,
                           legTitle=paste('P-value: ', pvalue),
                           title=nowotwor)$plot + 
          ylim(c(0,1)) + 
          xlim(c(0, max_time)) + 
          xlab("Time in days") + 
          ylab("Probability of survival") +
          theme(legend.position = c(0.9, 0.9))

  
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
      
    }, height = 800, width = 1000)
  
  output$opis_geny <- renderText({
    "The table shows the 10 most significant genes in which 
    a mutation occurred in cancer. \n"
  })
  
  output$geny <- renderTable({
    validate(
      need(input$nowotwory != "", "Please select a cancer!")
    )
    nowotwory <- input$nowotwory
   
      p <- matrix(1, nrow=297, ncol=1)
      for (nowotwor in nowotwory)
      {
        p2 <- p_value_tabela[order(p_value_tabela[,nowotwor]), ][, c("gen", nowotwor)]
        colnames(p2) <- c(paste(nowotwor,".Marker", sep=""), paste(nowotwor,".P-value", sep=""))
        p2[,2] <- signif(p2[,2], digits=5)
        p <- cbind(p, p2)
      }
      p <- p[, -c(1)]
      rownames(p) <- NULL
#       p <- ztable(p)
#       cgroup <- nowotwory
#       n.cgroup <- rep(2, length(nowotwory))
#       p <- addcgroup(p,cgroup=cgroup,n.cgroup=n.cgroup)
      print(p)
      
    }, digits = 5)
  
  output$opis_geny_wspol <- renderText({
    "The following table depicts three genes whose mutations appear most frequently with the mutation of a given gene in the given tumors."
  })
  
  output$geny_wspolne<-renderTable({
    validate(
      need(input$nowotwory != "", "Please select a cancer!")
    )
    gen <- input$geny
    nowotwory <- input$nowotwory
    
    #p <- matrix(1, nrow=3, ncol=1)
    p <- NULL
    for (nowotwor in nowotwory){
      dane <- get(paste('zbior.', nowotwor, sep=""))
      z <- numeric((ncol(dane)-2))
      for (i in 3:(ncol(dane)-2)){
        z[i] <- sum(dane[which(dane[, gen]==1), i])
        }
      a <- sum(dane[which(dane[, gen]==1), gen])
      x <- z[order(z)][(length(z)-1):(length(z)-3)]
      x2 <- c(which(z==x[1]), which(z==x[2]),which(z==x[3]))
      #x2 <- x2[1:3]
      
      x3 <- as.matrix(colnames(dane)[c(x2)])
      x4 <- as.matrix(round(z [ c(x2)]/a, 2))
      colnames(x3) <- c(paste("Marker ", nowotwor))
      colnames(x4) <- c(paste("Correlation ", nowotwor))
      p <- cbind(p, x3, x4)
      }
    #p <- p[, -c(1)]
    rownames(p) <- NULL
    print(p)
    })
  
  output$heatmap_pvalue <- renderPlot({
    gen <- input$geny
    #melted_dane <- melt(p_value_tabela[which(p_value_tabela$gen %in% geny[1:length(geny)]), ])
    melted_dane <- melt(p_value_tabela[which(p_value_tabela$gen==gen), ])
    base_size <- 12
    
    ggplot(data = melted_dane, aes(x=gen, y=variable, fill=value)) + 
      geom_tile() + theme_grey(base_size = base_size) + labs(x = "",y = "") + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(), 
                                                 axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))+
      geom_tile(aes(fill = value), colour = "white")
    
  })
  
  output$heatmap_czestosc <- renderPlot({
    melted_dane <- melt(czestosci[which(czestosci$gen %in% geny[1:length(geny)]), ])
    
    base_size <- 12
    
    
    ggplot(data = melted_dane, aes(x=variable, y=gen, fill=value)) + 
      geom_tile() + theme_grey(base_size = base_size) + labs(x = "",y = "") + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(), 
                                                 axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))+
      geom_tile(aes(fill = value), colour = "white") + 
      scale_fill_gradient(low = "white", high = "steelblue")
    
  })
  
  output$survcurves_variant <- renderPlot({
    validate(
      need(input$nowotwory != "", "Please select a cancer!")
    )
      nowotwory <- input$nowotwory
      gen <- input$geny
      
      n <- length(nowotwory)
      
      nowotwory_variant_all <- NULL
      for(nowotwor in nowotwory){
              
              
              nowotwor_variant <- get(paste(nowotwor, '_variant', sep=""))
              nowotwor_variant <- nowotwor_variant[, c("patient.barcode", "time", "status", paste('Variant.', gen, sep=""))]
              nowotwor_variant$nowotwor <- nowotwor
              colnames(nowotwor_variant) <- c("patient.barcode", "time", "status", "Variant", "nowotwor")
              #nowotwory.gen.missense <- rbind(nowotwory.gen.missense, nowotwor_variant[!is.na(nowotwor_variant$Variant) & nowotwor_variant$Variant== "Missense_Mutation", ])
              #nowotwory.gen.nonsense <- rbind(nowotwory.gen.nonsense, nowotwor_variant[!is.na(nowotwor_variant$Variant) & nowotwor_variant$Variant == "Nonsense_Mutation", ])
              
              nowotwor_variant$time <- as.numeric(as.character(nowotwor_variant$time))
              
              nowotwor_variant$missense <- as.numeric(ifelse(nowotwor_variant$Variant == "Missense_Mutation", 1, 0))
              nowotwor_variant$nonsense <- as.numeric(ifelse(nowotwor_variant$Variant == "Nonsense_Mutation", 1, 0))
              nowotwory_variant_all <- rbind(nowotwory_variant_all, nowotwor_variant[!is.na(nowotwor_variant$Variant),])
      }    
#       max_time <- 0
#       for(nowotwor in nowotwory){
#         zbior <- get(paste(nowotwor, '_variant', sep=""))
#         time <- as.numeric(as.character(zbior$time))
#         time <- max(time, na.rm = TRUE)
#         max_time <- ifelse(time>max_time, time, max_time)
#       } 
#       
      max_time <- max(nowotwory_variant_all$time, na.rm = TRUE)
      
      p.missense <- lapply(nowotwory, function(nowotwor){
        pvalue <- "Brak"
        
   
        dane <- nowotwory_variant_all[nowotwory_variant_all$nowotwor == nowotwor,]

        
 
      if(nrow(dane)>2){

        nowotwor_gen.fit.missense <- survfit(Surv(time, status) ~ missense, 
                                    data=dane)
        if (length(unique(dane$missense))==1){
          variant <- "No Missense"
        }
        else{
          variant <- c("No Missense", "Missense")
          
          survdiff <- survdiff(Surv(time, status) ~ missense, 
                               data=dane)
          pvalue <- signif(pchisq(survdiff$chisq, 1, lower=F), 3)
        }
        
        survMisc::autoplot(nowotwor_gen.fit.missense,
                           legLabs = variant,
                           legTitle=paste('P-value: ', pvalue),
                           title=paste(nowotwor, "\n  Missense Mutation", sep=""))$plot + 
          ylim(c(0,1)) + 
          xlim(c(0, max_time)) + 
          xlab("Time in days") + 
          ylab("Probability of survival") +
          theme(legend.position = c(0.85, 0.85))
      }
      })
      
        p.nonsense <- lapply(nowotwory, function(nowotwor){
        pvalue <- "Brak"
#         
#         nowotwory_variant_all <- NULL
#         nowotwor_variant <- get(paste(nowotwor, '_variant', sep=""))
#         nowotwor_variant <- nowotwor_variant[, c("patient.barcode", "time", "status", paste('Variant.', gen, sep=""))]
#         nowotwor_variant$nowotwor <- nowotwor
#         colnames(nowotwor_variant) <- c("patient.barcode", "time", "status", "Variant", "nowotwor")
#         #nowotwory.gen.missense <- rbind(nowotwory.gen.missense, nowotwor_variant[!is.na(nowotwor_variant$Variant) & nowotwor_variant$Variant== "Missense_Mutation", ])
#         #nowotwory.gen.nonsense <- rbind(nowotwory.gen.nonsense, nowotwor_variant[!is.na(nowotwor_variant$Variant) & nowotwor_variant$Variant == "Nonsense_Mutation", ])
#         nowotwory_variant_all <- rbind(nowotwory_variant_all, nowotwor_variant[!is.na(nowotwor_variant$Variant),])
#         nowotwory_variant_all$time <- as.numeric(as.character(nowotwory_variant_all$time))
#         
#         nowotwory_variant_all$missense <- as.numeric(ifelse(nowotwory_variant_all$Variant == "Missense_Mutation", 1, 0))
#         nowotwory_variant_all$nonsense <- as.numeric(ifelse(nowotwory_variant_all$Variant == "Nonsense_Mutation", 1, 0))
#         
        
        

        dane <- nowotwory_variant_all[nowotwory_variant_all$nowotwor == nowotwor,]

        if(nrow(dane)>2){

        nowotwor_gen.fit.nonsense <- survfit(Surv(time, status) ~ nonsense, 
                                             data=dane)
        
        if (length(unique(dane$nonsense))==1){
          variant <- "No Nonsense"
        }
        else{
          variant <- c("No Nonsense", "Nonsense")
          
          survdiff <- survdiff(Surv(time, status) ~ nonsense, 
                               data=dane)
          
          pvalue <- signif(pchisq(survdiff$chisq, 1, lower=F), 3)
        }
        
        survMisc::autoplot(nowotwor_gen.fit.nonsense,
                           legLabs = variant,
                           legTitle=paste('P-value: ', pvalue),
                           title=paste(nowotwor, "\n  Nonsense Mutation", sep=""))$plot + 
          ylim(c(0,1)) + 
          xlim(c(0, max_time)) + 
          xlab("Time in days") + 
          ylab("Probability of survival") +
          theme(legend.position = c(0.85, 0.85))
        }
      })
    
      for (i in 1:length(nowotwory))
      { 
        if (i==1)
        {
          z = p.missense[1]
        }
        else
        {
          z=append(z , p.missense[i] )
        }
        z = append(z, p.nonsense[i])
      }
      #marrangeGrob(append(p.missense, p.nonsense), nrow=length(nowotwory), ncol=2)
      marrangeGrob(z, nrow=length(nowotwory), ncol=2)
      
  }, height = 600, width = 1000)
  
  output$table_variant <- renderTable({
    validate(
      need(input$nowotwory != "", "Please select a cancer!")
    )
    nowotwor <- input$nowotwory
    gen <- input$geny
    
    p <- matrix(1, nrow=10, ncol=2)
    #p <- NULL
    for (nowotworr in nowotwor)
    {
      dane <- czestosci_variant[which(czestosci_variant$nowotwor==nowotworr), c("variant", gen)]
      colnames(dane) <- c(paste('variant ', nowotworr), paste("frequency in ", nowotworr))
      p <- cbind(p, dane)
    }
    p <- p[, -c(1, 2)]
    rownames(p) <- NULL
    print(p)
  }, digits=3)



  output$table_new <- renderDataTable({
  gen <- input$geny
  
  nowotwory_all <- c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                     "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
  dane<- matrix(0, nrow=14, ncol=4)
  p <- NULL


    dane[, 2] <- t(czestosci[czestosci$gen==gen,nowotwory_all])

    dane[,3]<-t(p_value_tabela[p_value_tabela$gen==gen, nowotwory_all])
  dane[, 1]<-nowotwory_all
    Licznosc<- c(1110, 1098, 941, 629, 628, 595, 591, 548, 537, 528, 522, 515, 504, 503)
    dane[,4]<-as.numeric(Licznosc)*as.numeric(dane[,2])
    colnames(dane)<-c('Cancer', 'Frequency', 'Pvalue', 'Cardinality of set')
 
  p<-cbind(p, dane)

  print(p)
}, options = list(dom = 't', lengthMenu = c(20, 30)))



  
  })