library(shiny)
library(survival)
library(survMisc)
#setwd<-'C:/Users/Gosia/Desktop/Biecek'
#setwd(setwd)
nowotwory <- list("Wszystkie", "GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

# dane = NULL
# for (nowotwor in nowotwory)
# {
#   if(nowotwor != "Wszystkie"){
#     #p_value\P_value_dla_interesujacych_genow
#     dane1 = read.table(paste('p_value/P_value_dla_interesujacych_genow/', nowotwor, '_pvalue.txt', sep=''), h=T)
#     dane1$gen = rownames(dane1)
#     dane1$nowotwor = nowotwor
#     dane = rbind(dane, dane1)
#   }
# }

#geny<-unique(dane$gen)

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)
#nowotwory <- list("Wszystkie", "GBMLGG", "BRCA", "KIPAN")

shinyUI(fluidPage(
  titlePanel("Mutacje genów"),
  sidebarLayout(
    sidebarPanel(
      selectInput("nowotwory",
                  "Wybierz nowotwór",
                  nowotwory,
                  "Wszystkie", multiple = TRUE),
      selectInput("geny",
                  "Wybierz gen",
                  geny,
                  "TP53")
    )
    ,
    
    mainPanel(
      p("Dane dotyczące wybranego mutacji wybranego genu:"),
      br(),
      tabsetPanel(
        tabPanel("Wykres", textOutput("text1"), plotOutput("wykres", width = 500)),
        tabPanel("Geny", tableOutput("geny")),
        tabPanel("Podsumowanie", verbatimTextOutput("podsumowanie"))        ,
        tabPanel("Tabela", tableOutput("tabela")),
        tabPanel("Kilka nowotworów", tableOutput("tabela_kilka"))
      )
    )
  )
)
)

