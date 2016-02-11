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
             For more information click'),
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
        tabPanel("Summary of gene mutation", 
                 HTML('<br/> For the selected gene
                 the following table contains information about the frequency and number 
                 of patients with mutation of this gene among patients suffering on the different types of cancers.
                 It also includes information about the importance of mutations on
                 patients survival measured by p-value of the log-rank test.<br/><br/>'),
                 dataTableOutput("table_new")),
        tabPanel("Survival curves: Presence of mutation", 
                 HTML('<br/>In the figures below we can see Kaplan-Meier curves for 
                 a given gene and the given tumors. The survival curves 
                 are estimated for the two groups of patients: 
                 the first one refers to the patients with a mutation of a given gene 
                 and the second one refers to the group of patients without any mutation 
                 of this gene.<br/><br/>'), plotOutput("survcurves_yesno")),
        tabPanel('Co-occuring genes', 
                 HTML('<br/> Opis <br/><br/>'), 
                 dataTableOutput("co_occuring_table")),
        tabPanel("Survival curves: Variant Classification", 
                 HTML('<br/> Opis <br/><br/>
                  <TABLE  WIDTH=100%>
                  <TR> 
                    <TD> <p align="center"> <font size="3"><b>Missense Mutation</b></font> </p> </TD> 
                    <TD> <p align="center"> <font size="3"><b>Nonsense Mutation</b></font> </p> </TD> 
                  </TR> </TABLE>
                      '), 
                 plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", 
                 HTML('<br/> Opis <br/><br/>'), 
                 dataTableOutput("table_variant"))
        )
      )
    )
  )
  )