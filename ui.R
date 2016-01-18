library(shiny)
library(survival)
library(survMisc)

nowotwory <- list("Wszystkie", "GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

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
    ),
    
    mainPanel(
      p("Dane dotyczące mutacji wybranego genu:"),
      br(),
      tabsetPanel(
        tabPanel("Wykres", textOutput("opis_wykres"), plotOutput("wykres", width = 500)),
        tabPanel('Geny współwystępujące', textOutput("opis_geny_wspol"), tableOutput("geny_wspolne")),
        tabPanel("Geny", textOutput("opis_geny"), tableOutput("geny"))
      )
    )
  )
)
)