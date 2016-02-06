library(shiny)
library(survival)
# library(survMisc)
library(stats)

nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

shinyUI(fluidPage(
  titlePanel("Genes mutations and survival analysis"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("nowotwory",
                  "Select cancer",
                  nowotwory,
                  "BRCA", multiple = TRUE, options = list(maxItems = 4)),
      selectInput("geny",
                  "Select gene",
                  geny,
                  "TP53")
    ),
    
    mainPanel(
      p(""),
      br(),
      tabsetPanel(
        tabPanel("Survival curves: presence of mutation", textOutput("opis_krzywe"), plotOutput("survcurves_yesno")),
        tabPanel('Co-occuring genes', textOutput("opis_geny_wspol"), tableOutput("geny_wspolne")),
        tabPanel("Most significant genes", textOutput("opis_geny"), tableOutput("geny")),
        tabPanel("Heatmap: p-value", plotOutput("heatmap_pvalue", width = 600, height = 2500)),
        tabPanel("Heatmap: frequency", plotOutput("heatmap_czestosc", width = 600, height = 2500)),
        tabPanel("Survival curves: Variant Classification", plotOutput("survcurves_variant", width = 500, height = 700)),
        tabPanel("Table: Variant Classification", tableOutput("table_variant")),
        tabPanel("nowa", dataTableOutput("table_new"))
        
      )
    )
  )
)
)