library(shiny)
library(survival)
# library(survMisc)
library(stats)
library(DT)

nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

shinyUI(fluidPage(
  titlePanel(div(img(src="logo.png", height = 60, width = 60), "Genes mutations")),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("nowotwory",
                  "Select cancer",
                  nowotwory,
                  "BRCA", multiple = TRUE, options = list(maxItems = 4)),
      selectInput("geny",
                  "Select gene",
                  geny,
                  "TP53"),
      tags$div(
        HTML('<br/><br/><br/>
             <font size="2"><b>Details:</b><br/>
             <br/><br/>
             <b>Authors:</b><br/>
             Sebastian Gargas, Małgorzata Dobkowska, Marlena Bielat</font>')
      )
    ),
    
    mainPanel(
      p(""),
      br(),
      tabsetPanel(
        tabPanel("nowa", textOutput('opis_nowa'),dataTableOutput("table_new")),
        tabPanel("Survival curves: Presence of mutation", textOutput("opis_krzywe"), plotOutput("survcurves_yesno")),
        tabPanel("Survival curves: Variant Classification", plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", dataTableOutput("table_variant")),
        tabPanel('Co-occuring genes', textOutput("opis_geny_wspol"),
                 hr(),
                 dataTableOutput("geny_wspolne"))
        #tabPanel("Most significant genes", textOutput("opis_geny"), tableOutput("geny")),
        #tabPanel("Heatmap: p-value", plotOutput("heatmap_pvalue", width = 600, height = 2500)),
        #tabPanel("Heatmap: frequency", plotOutput("heatmap_czestosc", width = 600, height = 2500)),
      )
    )
  )
)
)