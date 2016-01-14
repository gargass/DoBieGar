library(shiny)
library(survival)
library(survMisc)
#setwd<-'C:/Users/Gosia/Desktop/Biecek'
#setwd(setwd)
nowotwory <- list("Wszystkie", "GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

dane = NULL
for (nowotwor in nowotwory)
{
  if(nowotwor != "Wszystkie"){
    dane1 = read.table(paste('pvalue/', nowotwor, '_pvalue.txt', sep=''), h=T)
    dane1$gen = rownames(dane1)
    dane1$nowotwor = nowotwor
    dane = rbind(dane, dane1)
  }
}
geny<-unique(dane$gen)


#nowotwory <- list("Wszystkie", "GBMLGG", "BRCA", "KIPAN")

shinyUI(fluidPage(
  titlePanel("Mutacje genów"),
  sidebarLayout(
    sidebarPanel(
      selectInput("geny",
                  "Wybierz gen",
                  geny,
                  "TP53"),
      selectInput("nowotwory",
                  "Wybierz nowotwór",
                  nowotwory,
                  "Wszystkie")
    )
    ,
    
    mainPanel(
      p("Dane dotyczące wybranego mutacji wybranego genu:"),
      br(),
      tabsetPanel(
        tabPanel("Wykres", plotOutput("wykres", width = 500)),
        tabPanel("Podsumowanie", verbatimTextOutput("podsumowanie"))        ,
        tabPanel("Tabela", tableOutput("tabela")),
        tabPanel("Kilka nowotworów", tableOutput("tabela_kilka"))
      )
    )
  )
)
)

