library(shiny)
library(survival)
# library(survMisc)
library(stats)
library(DT)

nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

names(nowotwory) <- c("GBMLGG - Glioblastoma multiforme",
                    "BRCA - Breast invasive carcinoma",
                    "KIPAN - Pan-kidney cohort",
                    "COADREAD - Colorectal adenocarcinoma",
                    "STES - Stomach and Esophageal carcinoma",
                    "GBM - Glioblastoma multiforme",
                    "OV - Ovarian serous cystadenocarcinoma",
                    "UCEC - Uterine Corpus Endometrial Carcinoma",
                    "KIRC - Kidney renal clear cell carcinoma",
                    "HNSC - Head-Neck Squamous Cell Carcinoma",
                    "LUAD - Lung adenocarcinoma",
                    "LGG - Lower Grade Glioma",
                    "LUSC - Lung squamous cell carcinoma",
                    "THCA - Thyroid carcinoma")

geny <- read.table('p_value/lista_interesujacych_genow.txt', h=T)
geny <- as.matrix(geny)

shinyUI(fluidPage(
  titlePanel(div(img(src="logo.png", height = 100, width = 100), "Genes mutations")),
  sidebarLayout(
    sidebarPanel(
      selectInput("geny",
                  "Select gene",
                  geny,
                  "TP53"),
      selectizeInput("nowotwory",
                  "Select cancer",
                  nowotwory,
                  "BRCA", multiple = TRUE, options = list(maxItems = 4)),
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
      tabsetPanel(
        tabPanel("Instruction", 
                 tags$h2("Instrukcja"),
                 HTML(' Hello World <br/> <b>hahaha</b> <br/>
                      <font size="5"> halo halo </font>
                      <ul> 
                        <li>czesc</li> 
                        <li> hi </li> 
                     </ul>
                      ')),
        tabPanel("Basic information about the gene mutation", textOutput('basic_description'),dataTableOutput("table_new")),
        tabPanel("Survival curves: Presence of mutation", textOutput("curves_description"), plotOutput("survcurves_yesno")),
        tabPanel('Co-occuring genes', textOutput("co_occuring_description"),hr(), dataTableOutput("co_occuring_table")),
        tabPanel("Survival curves: Variant Classification", plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", dataTableOutput("table_variant"))

        )
      )
    )
  )
  )