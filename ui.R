library(shiny)
library(survival)
library(survMisc)

nowotwory <- list("Wszystkie", "GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

shinyUI(fluidPage(
  titlePanel("Gene Mutation"),
  sidebarLayout(
    sidebarPanel(
      selectInput("nowotwory",
                  "Select cancer",
                  nowotwory,
                  "Wszystkie", multiple = TRUE),
      selectInput("geny",
                  "Select gene",
                  geny,
                  "TP53")
    ),
    
    mainPanel(
      p("Dane dotyczące mutacji wybranego genu:"),
      br(),
      tabsetPanel(
        tabPanel("Survival curve", textOutput("opis_krzywe"), plotOutput("wykres", width = 500)),
        tabPanel('Co-occuring genes', textOutput("opis_geny_wspol"), tableOutput("geny_wspolne")),
        tabPanel("10 most significant genes", textOutput("opis_geny"), tableOutput("geny"))
      )
    )
  )
)
)