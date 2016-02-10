library(shiny)
library(survival)
# library(survMisc)
library(stats)
library(DT)

nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

names(nowotwory) <- c("GBMLGG - Glioblastoma multiforme",
                    "BRCA - Breast invasive carcinoma",
                    "KIPAN",
                    "COADREAD - Colorectal adenocarcinoma",
                    "STES",
                    "GBM",
                    "OV - Ovarian serous cystadenocarcinoma",
                    "UCEC - Uterine Corpus Endometrial Carcinoma",
                    "KIRC - Kidney renal clear cell carcinoma",
                    "HNSC",
                    "LUAD - Lung adenocarcinoma",
                    "LGG - Lower Grade Glioma",
                    "LUSC - Lung squamous cell carcinoma",
                    "THCA")

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
             For more information click '),
        tags$a("here",target="_blank",href="doc.pdf"),
        HTML('.<br/><br/>
             <b>Authors:</b><br/>
             <a href="mailto:gargass@student.mini.pw.edu.pl">Sebastian Gargas</a><br/>
             <a href="mailto:dobkowskam@student.mini.pw.edu.pl">Małgorzata Dobkowska</a><br/>
             <a href="mailto:bielatm@student.mini.pw.edu.pl">Marlena Bielat</a><br/></font>')
      )
    ),
    
    mainPanel(
      #p(""),
      #br(),
      tabsetPanel(
        tabPanel("nowa", textOutput('opis_nowa'),dataTableOutput("table_new")),
        tabPanel("Survival curves: Presence of mutation", textOutput("opis_krzywe"), plotOutput("survcurves_yesno")),
        tabPanel("Survival curves: Variant Classification", plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", dataTableOutput("table_variant")),
        tabPanel('Co-occuring genes', textOutput("opis_geny_wspol"),
                 hr(),
                 dataTableOutput("geny_wspolne"))
        )
      )
    )
  )
  )