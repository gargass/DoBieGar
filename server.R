#Loading packages
library(shiny)
library(survival)
# library(survMisc)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(DT)
library(scales)

#Loading datasets
nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

for(nowotwor in nowotwory){
  assign(paste('zbior.', nowotwor, sep=""), read.table(paste('Zbiory/', nowotwor, '.txt', sep="")))
  assign(paste('geny_wspolne_', nowotwor, sep=""), read.table(paste('Zbiory/wspolne_', nowotwor, '.txt', sep="")))
  assign(paste('geny_wspolne_licznosci_', nowotwor, sep=""), read.table(paste('Zbiory/wspolne_', nowotwor, '_licznosci.txt', sep="")))
  
}

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

p_value_tabela <-read.table('p_value/p_value_NA.txt', h=T)
czestosci<-read.table('czestosci.txt', h=T)
licznosci<-read.table('licznosci.txt', h=T)

for(nowotwor in nowotwory){
  assign(paste('p_value.', nowotwor, sep=""), read.table(paste('p_value/P_value_dla_interesujacych_genow/', 
                                                               nowotwor, '_pvalue.txt', sep=""), h=T))
}

for(nowotwor in nowotwory){
  assign(paste(nowotwor, '_variant', sep=""), read.table(paste('Zbiory/', nowotwor, '_variant.txt', sep="")))
}

czestosci_variant <- read.table("czestosci_variant.txt", h=T)



#Shiny
shinyServer(function(input, output) {

#Basic information about gene mutation

#Table
  output$table_new <- renderDataTable({
    gen <- input$geny
    
    nowotwory_all <- c("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                       "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")
    
    
    freq <- round(as.numeric(czestosci[czestosci$gen==gen,nowotwory_all]),3)
    
    dane<- data.frame(matrix(0, nrow=14, ncol=4))
    colnames(dane) <- c("cancer", "freq", "n", "pvalue")
    dane$cancer <-nowotwory_all
    #dane$freq <- percent(freq)
    dane$freq<-freq
    dane$n <- as.numeric(licznosci[licznosci$gen==gen, nowotwory_all])
    pv<-as.numeric(p_value_tabela[p_value_tabela$gen==gen, nowotwory_all])
    for (i in 1:length(nowotwory))
    {
      pv[i]<-signif(pv[i], 3)
    }
    dane$pvalue <- pv
#       dane <- matrix(1, nrow=14, ncol=4)
  
#       dane[,1]<-t(nowotwory_all)
#       dane[,2]<-t(as.numeric(round(as.numeric(czestosci[czestosci$gen==gen,nowotwory_all]),3)))
#       dane[,3]<-as.numeric(t(licznosci[licznosci$gen==gen, nowotwory_all]))
#       dane[,4]<-t(as.numeric(round(p_value_tabela[p_value_tabela$gen==gen, nowotwory_all], 6)))
# 
#     
    colnames(dane)<-c('Cancer', 'Mutation frequency', 'Number of patients with mutation', 'Significance')
    dane
  }, options = list(columnDefs= list(list(className = 'dt-right', targets='_all')), dom = 't', lengthMenu = c(20, 30)))


#Survival Curves - Presence of mutation
  
#Curves
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
    
    survMisc::autoplot(nowotwor_gen.fit, legLabs = mutation,
                       legTitle=paste('P-value: ', pvalue),
                       title=nowotwor, censSize=2)$plot + 
      ylim(c(0,1)) + 
      xlim(c(0, max_time)) + 
      xlab("Time in days") + 
      ylab("Survival") +
      theme(legend.position = c(0.85, 0.9)) + 
      theme(legend.title = element_text(colour=ifelse(pvalue<0.05,"red", "black"), face="bold"))
  })
  
  if(n <= 4){
    ncol <- 2
    nrow <- 2
  }
  else{
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n/ceiling(sqrt(n)))
  }
  
  marrangeGrob(p, ncol = ncol, nrow = nrow, top=NULL)
  
}, height = 800)


#Co-occuring genes




#Table

output$co_occuring_table<-renderDataTable({
  validate(
    need(input$nowotwory != "", "Please select a cancer!")
  )
  gen <- input$geny
  nowotwory <- input$nowotwory
  
  tabela<-NULL
  for(nowotwor in nowotwory){
    
    gen_x_gen <- get(paste('geny_wspolne_', nowotwor, sep=""))
    geny <- rownames(gen_x_gen)[-which(rownames(gen_x_gen) == gen)]
    freq <- round(as.numeric(gen_x_gen[geny, gen]),2)
    #tabela <- cbind(tabela, freq)
    gen_y_gen<-get(paste('geny_wspolne_licznosci_', nowotwor, sep=""))
    #tabela<-cbind(tabela, gen_y_gen[geny, gen])
    n <- gen_y_gen[geny, gen]
    freq_n <- paste(freq, ' (', n, ')', sep="")
    tabela <- cbind(tabela, freq_n)
 
  }

  rownames(tabela) <- geny
  col<-NULL
  for (i in 1:length(nowotwory))
  {
  col <- append(col,nowotwory[i])
  #col <- append(col,paste(nowotwory[i], 'Number of patiens'))
  }
  colnames(tabela)<-col
  tabela
}, options = list( columnDefs = list(list(className = 'dt-right', targets='_all'))),filter='bottom')



#Variant Classification
#Curves
  
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
              
        nowotwor_variant$time <- as.numeric(as.character(nowotwor_variant$time))
              
        nowotwor_variant$missense <- as.numeric(ifelse(nowotwor_variant$Variant == "Missense_Mutation", 1, 0))
        nowotwor_variant$nonsense <- as.numeric(ifelse(nowotwor_variant$Variant == "Nonsense_Mutation", 1, 0))
        nowotwory_variant_all <- rbind(nowotwory_variant_all, nowotwor_variant[!is.na(nowotwor_variant$Variant),])
      }    
 
      max_time <- max(nowotwory_variant_all$time, na.rm = TRUE)
      
      p.missense <- lapply(nowotwory, function(nowotwor){
        pvalue <- "Brak"
        dane <- nowotwory_variant_all[nowotwory_variant_all$nowotwor == nowotwor,]
        
        if(nrow(dane)>2){
          nowotwor_gen.fit.missense <- survfit(Surv(time, status) ~ missense, data=dane)
          if (length(unique(dane$missense))==1){
            variant <- "No Missense"
            }
          else{
            variant <- c("No Missense", "Missense")
            survdiff <- survdiff(Surv(time, status) ~ missense, data=dane)
            pvalue <- signif(pchisq(survdiff$chisq, 1, lower=F), 3)
            }
          
          survMisc::autoplot(nowotwor_gen.fit.missense,
                           legLabs = variant,
                           legTitle=paste('P-value: ', pvalue),
                           #title=paste(nowotwor, "\n  Missense Mutation", sep=""), censSize=2)$plot +
                           title=nowotwor, censSize=2)$plot + 
            ylim(c(0,1)) + 
            xlim(c(0, max_time)) + 
            xlab("Time in days") + 
            ylab("Survival") + 
            theme(legend.position = c(0.85, 0.9)) + 
            theme(legend.title = element_text(colour=ifelse(pvalue<0.05,"red", "black"), face="bold"))
          }
        })
      
      p.nonsense <- lapply(nowotwory, function(nowotwor){
        pvalue <- "NULL"

        dane <- nowotwory_variant_all[nowotwory_variant_all$nowotwor == nowotwor,]

        if(nrow(dane)>2){
          nowotwor_gen.fit.nonsense <- survfit(Surv(time, status) ~ nonsense, data=dane)
          
          if (length(unique(dane$nonsense))==1){
            variant <- "No Nonsense"
            }
          else{
            variant <- c("No Nonsense", "Nonsense")
            survdiff <- survdiff(Surv(time, status) ~ nonsense, data=dane)
            pvalue <- signif(pchisq(survdiff$chisq, 1, lower=F), 3)
            }
          
          survMisc::autoplot(nowotwor_gen.fit.nonsense,
                           legLabs = variant,
                           legTitle=paste('P-value: ', pvalue),
                           #title=paste(nowotwor, "\n  Nonsense Mutation", sep=""), censSize=2)$plot + 
                           title=nowotwor, censSize=2)$plot + 
            ylim(c(0,1)) + 
            xlim(c(0, max_time)) + 
            xlab("Time in days") + 
            ylab("Survival") + 
            theme(legend.position = c(0.85, 0.9)) + 
            theme(legend.title = element_text(colour=ifelse(pvalue<0.05,"red", "black"), face="bold"))
          }
        })
      
      for (i in 1:length(nowotwory)){ 
        if (i==1){
          z = p.missense[1]
        }
        else{
          z=append(z , p.missense[i] )
        }
        z = append(z, p.nonsense[i])
      }
      indeks <- NULL
      k <- 1
      for(e in z){
        if(is.null(e[[1]])){indeks <- append(indeks, k)
        }
        k <- k+1
        }
      if(!is.null(indeks)){
        z <- z[-indeks]
        }
      if(length(z)>0){
        marrangeGrob(z, nrow=length(z)/2, ncol=2, top=NULL)
      }
      else{
        validate(
          need(length(z)>0, "Too few observations in the selected group!")
        )
      }
      }, height = 800)
  

#Variant Classification
#Table
  output$table_variant <- renderDataTable({
    validate(
      need(input$nowotwory != "", "Please select a cancer!")
    )
    nowotwor <- input$nowotwory
    gen <- input$geny

    p <- matrix(1, nrow=12, ncol=length(nowotwor))
    k <- 1

    if(gen %in% colnames(czestosci_variant)){
      for (nowotworr in nowotwor){
        dane <- czestosci_variant[which(czestosci_variant$nowotwor==nowotworr), c("variant", gen)]
        p[3:12, k] <- dane[,2]
        
        k <- k + 1
        }
    }
    else{
      for (nowotworr in nowotwor){
        p[3:12, k] <- rep(0, 10)
        
        k <- k + 1
      }
    }
    
    for (k in 1:(length(nowotwor)))
    {
      if (sum(as.numeric(p[3:12,k]))!=0)
      {
   
      p[3:12,k]= paste(round((100*as.numeric(p[3:12,k]))/sum(as.numeric(p[3:12,k])),2), "%", sep="")
      }
      else
      {
        p[3:12,k]= rep(paste(0, "%", sep=""), 10)
      }
    }
    #dane[, 2] <- t(paste(round(100*czestosci[czestosci$gen==gen,nowotwory_all],3), "%", sep=""))
    
    p[1,1:(length(nowotwor))]<-t(licznosci[licznosci$gen==gen, nowotwor])
    p[2,1:(length(nowotwor))] <- paste(round(100*czestosci[czestosci$gen==gen,nowotwor],3), "%", sep="")
    
    rownames(p) <- c('Number of patients with mutation', 'Mutation frequency','Missense_Mutation', 'Silent', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'Nonsense_Mutation', 
               'RNA', 'Splice_Site', 'In_Frame_Ins', 'Nonstop_Mutation')

    colnames(p) <- c(nowotwor)
    p
  },  options = list(columnDefs = list(list(className = 'dt-right', orderable= FALSE)), dom='t', paging=FALSE))



  })