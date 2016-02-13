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
             Full documentation is available '),
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
                        <li> 14 types of cancer with data from at least 500 patients</li> 
                        <li> 535 biomarkers - mutations on specific genes. </li> 
                        <li> Each gene considered by us has significant effect on 
                             survival of patients suffering from at least one of considered cancers. </li>
                        <li> Keep in mind that mutations of genes are rare events.</li>
                     </ul>

                      <br/>
                    
                      '),
                 
                 tags$h4("What you can find here?"),
                 HTML('
                      At the beginning, select one of biomarkers which interests you, then select cancers for which you 
                      want to check gene mutation impact on survival time. You can choose at most 4 cancer types.
                      On the another panels you can find following informations:
                      <ul>
                        <li> Summary of gene mutation : 
                          <ul>
                                <li> Frequency of mutation for each cancers</li>
                                <li> Number of patients with mutation for each cancers</li>
                                <li> Significance of mutations to the patients survival for each cancers</li>
                          </ul>
                        </li>
                      <li> Survival curves: Presence of mutation: </li>
                          <ul>
                                <li> Kaplan-Meier curves estimated for two groups of patients: 
                                with and without any mutation of the selected gene in selected cancers.</li>
                           </ul>
                      </li>
                      <li>Co-occuring genes
                          <ul>
                            <li>List of all considered biomarkers with co-occurence with mutation on the selected gene measured by
                          frequency patients with mutation on biomarker in row among patients with mutation on the selected gene
                          and number of patients with both mutations in specified cancer types.</li>
                          </ul>
                      </li>
                      <li>Survival curves: Variant Classification
                        <ul>
                          <li>Kaplan-Meier curves estimated for the two 
                          groups of patients: with a specific type of mutation - Missense Mutation or Nonsense Mutation 
                          and with the other type of mutation of this gene.
                        </ul>
                      </li>
                      <li>Frequency of mutation types
                        <ul>
                          <li>Occurence of frequency of different types of mutations for 
                              a selected gene and selected cancers.</li>
                        </ul>
                      </li>
                    </ul>
                      '),
                 
  
                 
                 HTML('###########################################<br/>'),
                 
                 
                 
                 HTML('It was created during the course of the subject 
                      Research Workshop at the Faculty of Mathematics and 
                      Information Science at Warsaw University of Technology. 
                      Its purpose is to show the effect of the mutation of 
                      different genes on survival of patients with various 
                      types of cancer. All analyzes are based on data from T
                      he Cancer Genome Atlas. Mutations are rare events, so it 
                      is worth keep in mind when you watch the results.<br/><br/>
                      <font size="4">How to start?</font><br/>
                      At the beginning choose the gene that interests you. 
                      In the bookmark "Summary of gene mutation" you find 
                      information about the frequency and the number of 
                      patients who have had a mutation in the selected gene 
                      among patients with 14 types of cancer, which we consider 
                      in our application. Here you can also see whether a mutation 
                      in selected gene significantly affects the survival time 
                      of patients. This significance is measured by p-value of 
                      the log-rank test comparing survival among patients who 
                      have had the mutation and those who did not have the 
                      mutation.<br/><br/>
                      In the bookmark "Survival Curves - Presence of the mutation" 
                      you also have to choose cancer for which you want to see 
                      an estimate of the probability of survival using the 
                      Kaplan Meier method. You can choose up to 4 cancer. 
                      K-M curves are estimated in two groups of patients with 
                      a particular type of cancer - in which has occurred and 
                      there was no mutation of the gene. It is worth noting that 
                      the presence of the mutation may have significant positive 
                      or negative effects on the survival time of patients: 
                      e.g., for the presence of the TP53 gene, the mutation 
                      in this gene has a positive effect on cancer GBMLGG 
                      and negatively on cancer STES.<br/><br/>'),
                 HTML(' Hello World <br/> <b>hahaha</b> <br/>
                      <font size="5"> halo halo </font>
                      <ul> 
                        <li>czesc</li> 
                        <li> hi </li> 
                     </ul>
                      '),
                HTML('Na razie po polsku<br/>
                     Witamy w naszej aplikacji! Powstała ona podczas zajęc z przedmiot Wasztaty Badawcze
                     prowadzonych na Wydziale Matematyki i Nauk Informacyjnych Politechniki Warszawskiej.
                     Jej celem jest pokazanie wpływu mutacji różnych genów na czas przeżycia pacjentów 
                     cierpiących na różne typy nowotworów. Wszystkie analizy oparte są na danych pochodzących 
                     z Genome Atlas tralala. Mutacje są rzadko występującymi zdarzeniami, więc warto podczas 
                     oglądania rezultatów mieć to na uwadze.

                    <br/>

                    Jak rozpocząć?
                    <br/>
                    Na początku wybierz interesujący Cię gen.
                    W zakładce Summary of gene mutation znajdziesz informację o częstości 
                    i liczbie pacjentów u których wystąpiła mutacja na wybranym
                    genie wśród pacjentów cierpiących na 14 typów nowotworów, które rozważamy w naszej aplikacji.
                    Również możecie tu zobaczyć, czy mutacja na wybranych genie istotnie wpływa na czas przeżycia 
pacjentów. Istotność ta mierzona jest p-wartością testu log-rank porównującego czas przeżycia wśród pacjentów u których
wystąpiła mutacja i tych, którzy danej mutacji nie mieli.
                    <br/>

W zakładce Survival Curves - Presence of mutation należy oprócz genu wybrać również nowotwory, dla których chcemy
zobaczyc oszacowanie prawdopodobieństwa przeżycia za pomocą krzywej Kaplana Meiera. Możesz wybrac maksymalnie 4 nowotwory.
Krzywe K-M są szacowane w dwóch grupach pacjentów cierpiącyh na dany typ nowotworu - u których wystapiła i niewystąpiła 
mutacja danego genu. Warto zwrócić uwagę, że obecność mutacji może mieć isttny pozytywny i negatywny wplyw na czas przezycia 
pacjentow:
np dla genu TP53 obncność ,mutacji na tym genie wpływa pozytywnie na chorych na raka GBMLGG, podczas gdy na pacjeentow 
cierpiacyhc na SYES wypływa negatywnie.

                     ')),
        
        
        
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
                 HTML('<br/> 

The table below presents informations about co-occurence mutations. 
For the selected gene and selected cancers a list of all considered biomarkers is given  
with information about the occurrence of mutations of these genes among
patients with mutation the selected gene in selected cancers types. In parentheses are given the cardinality of patients with the mutation in selected genes.
<br/><br/>'), 
                 dataTableOutput("co_occuring_table")),
        
        
        tabPanel("Survival curves: Variant Classification", 
                 HTML('<br/> 
The following figures present the Kaplan-Meier curves 
for the given gene and given cancers, for the two main types 
of mutation. The survival curves are estimated for the two 
groups of patients: the first one refers to the patients with a specific 
type of mutation - Missense Mutation or Nonsense Mutation and the second 
one refers to the patients with the other type of mutation of this gene.

<br/><br/>
                  <TABLE  WIDTH=100%>
                  <TR> 
                    <TD> <p align="center"> <font size="3"><b>Missense Mutation</b></font> </p> </TD> 
                    <TD> <p align="center"> <font size="3"><b>Nonsense Mutation</b></font> </p> </TD> 
                  </TR> </TABLE>
                      '), 
                 plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", 
                 HTML('<br/> 
The table below contains informations
about the occurrence frequency 
of different types of mutations for a selected gene and selected cancers.
In addition, we again show the level of occurrence of mutations in a given gene.
                      <br/><br/>'), 
                 dataTableOutput("table_variant"))
        )
      )
    )
  )
  )