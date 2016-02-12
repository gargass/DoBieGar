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
                 tags$h2("Welcome to our application!"),
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
                 The following table contains information about the frequency and number 
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
                 HTML('<br/> The table below presents information about co-occurence mutations. For selected gene
                      and selected cancers we can see list of all consider genes and information
                      about occurence of mutations of these genes among patients
                      with mutations of selected gene in selected cancers types (number of patients in parentheses).<br/><br/>'), 
                 dataTableOutput("co_occuring_table")),
        
        
        tabPanel("Survival curves: Variant Classification", 
                 HTML('<br/> In the figures below we can see Kaplan-Meier curves for 
                 a given gene and the given tumors for two main types of mutation. The survival curves 
                 are estimated for the two groups of patients: 
                      the first one refers to the patients with a specific type of mutation: Missense Mutation
                      or Nonsense Mutation
                      and the second one refers to the patients with other type of mutation 
                      of this gene. <br/><br/>
                  <TABLE  WIDTH=100%>
                  <TR> 
                    <TD> <p align="center"> <font size="3"><b>Missense Mutation</b></font> </p> </TD> 
                    <TD> <p align="center"> <font size="3"><b>Nonsense Mutation</b></font> </p> </TD> 
                  </TR> </TABLE>
                      '), 
                 plotOutput("survcurves_variant")),
        tabPanel("Frequency of mutation types", 
                 HTML('<br/> For the selected gene and selected cancers the table below contains information
                      about occurance frequency of different types of mutations. In addition, we show again 
                      the level of occurence of mutations in a given gene.
                      <br/><br/>'), 
                 dataTableOutput("table_variant"))
        )
      )
    )
  )
  )