library(shiny)
library(survival)
# library(survMisc)
library(stats)
library(DT)

nowotwory <- list("GBMLGG", "BRCA", "KIPAN", "COADREAD", "STES", "GBM", "OV",
                  "UCEC", "KIRC", "HNSC", "LUAD", "LGG", "LUSC", "THCA")

names(nowotwory) <- c("GBMLGG - GBM + LGG",
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
             Full documentation is available '),
        tags$a("here",target="_blank",href="doc.pdf"),
        HTML('.<br/><br/>
             <b>Authors:</b><br/>
             <a href="mailto:bielat.marlena@gmail.com">Marlena Bielat</a><br/>
             <a href="mailto:gosia.dobkowska@tlen.pl">Małgorzata Dobkowska</a><br/>
             <a href="mailto:gargass@student.mini.pw.edu.pl">Sebastian Gargas</a><br/></font>')
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Instruction", 
                 tags$h3("Welcome!"),
                 
                 HTML('This aplication was created for research workshop conducted at the 
                       Faculty of Mathematics and Information Science at Warsaw University of Technology. 
                       This aplication allows you to browse the results of analysis of mutations of genes
                       from The Cancer Genome Atlas 
                       (<a href="http://cancergenome.nih.gov">http://cancergenome.nih.gov</a>) 
                       using RTCGA package 
                       (<a href="https://github.com/RTCGA/RTCGA.mutations">https://github.com/RTCGA/RTCGA.mutations</a>).
                      <br/><br/>
                      '),
                 
                 tags$h4("Specification of data"),
                 HTML('<br/>
                     <ul> 
                        <li> 14 types of cancer with data from at least 500 patients.</li> 
                        <li> 535 biomarkers - mutations on specific genes. </li> 
                        <li> Each gene considered by us has significant effect on 
                             survival of patients suffering from at least one of considered cancers. </li>
                        <li> Keep in mind that mutations of genes are rare events.</li>
                     </ul><br/>'),
                 
                 tags$h4("What you can find here?"),
                 HTML('
                      At the beginning, select one of biomarkers which interests you, then select cancers for which you 
                      want to check gene mutation impact on survival time. You can choose at most 4 cancer types.
                      On the another panels you can find following informations:
                      <ul>
                        <li> Summary of gene mutation: 
                          <ul>
                                <li> Frequency of mutation for each cancers.</li>
                                <li> Number of patients with mutation for each cancers.</li>
                                <li> Significance of mutations to the patients survival for each cancers.</li>
                          </ul>
                        </li>
                      <li> Survival curves: Presence of mutation: </li>
                          <ul>
                                <li> Kaplan-Meier curves estimated for two groups of patients: 
                                with and without any mutation of the selected gene in selected cancers.</li>
                           </ul>
                      </li>
                      <li>Co-occuring genes:
                          <ul>
                            <li>List of all considered biomarkers with co-occurence with mutation on the selected gene measured by
                          frequency patients with mutation on biomarker in row among patients with mutation on the selected gene
                          and number of patients with both mutations in specified cancer types.</li>
                          </ul>
                      </li>
                      <li>Survival curves: Variant Classification:
                        <ul>
                          <li>Kaplan-Meier curves estimated for the two 
                          groups of patients: with a specific type of mutation - Missense Mutation or Nonsense Mutation 
                          and with the other type of mutation of this gene.
                        </ul>
                      </li>
                      <li>Frequency of mutation types:
                        <ul>
                          <li>Occurence of frequency of different types of mutations for 
                              a selected gene and selected cancers.</li>
                        </ul>
                      </li>
                    </ul>')),
        
        tabPanel("Summary of gene mutation", 
                 HTML('<br/> 
                 The following table contains informations about the frequency and number 
                 of patients with the mutation of the selected gene among patients suffering from different types of cancers.
                 It also includes information about the significance of mutations to the
                 patients survival measured by the p-value of the log-rank test.<br/><br/>'),
                 dataTableOutput("table_new")),
        tabPanel("Survival curves: Presence of mutation", 
                 HTML('<br/>The figures below ilustrate Kaplan-Meier curves for 
                 the given gene and the given cancers. The survival curves 
                 are estimated for the two groups of patients: 
                 the first one refers to the patients with a mutation of the given gene 
                 and the second one refers to the group of patients without any mutation 
                 of this gene.<br/><br/>'), plotOutput("survcurves_yesno")),
        tabPanel('Co-occuring genes', 
                 HTML('<br/> The table below presents informations about 
                  co-occurence mutations. For the selected gene and selected 
                  cancers a list of all considered biomarkers is given  
                  with information about the occurrence of mutations of these 
                  genes among patients with mutation the selected gene in selected 
                  cancers types. In parentheses are given the number of 
                  patients with the mutation in selected genes.<br/><br/>'), 
                 dataTableOutput("co_occuring_table")),
        
        
        tabPanel("Survival curves: Variant Classification", 
                 HTML('<br/> The following figures present the Kaplan-Meier curves 
                    for the given gene and given cancers, for the two main types 
                    of mutation. The survival curves are estimated for the two 
                    groups of patients: the first one refers to the patients with 
                    a specific type of mutation - Missense Mutation or Nonsense 
                    Mutation and the second one refers to the patients with the 
                    other type of mutation of this gene.
                  <br/><br/>
                  <TABLE  WIDTH=100%>
                  <TR> 
                    <TD> <p align="center"> <font size="3"><b>Missense Mutation</b></font> </p> </TD> 
                    <TD> <p align="center"> <font size="3"><b>Nonsense Mutation</b></font> </p> </TD> 
                  </TR> </TABLE>'), 
                 plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", 
                 HTML('<br/> The table below contains informations
                      about the occurrence frequency of different types of 
                      mutations for a selected gene and selected cancers.
                      In addition, we again show the level of occurrence of 
                      mutations in a given gene.<br/><br/>'), 
                 dataTableOutput("table_variant"))
        )
      )
    )
  )
  )