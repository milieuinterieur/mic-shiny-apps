
# 
# Milieu Interieur Consortium
# Companion Shiny App Nanostring25 
# 
# Author: Vincent Rouilly
# Oct 2015
#
# UI.R 

library(shiny)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("Standardized whole-blood transcriptional profiling enables the deconvolution of complex induced immune responses"),
    
    h4("Interactive Application provided by", align = "center"),
    h4(a("The Milieu Interieur Consortium", href="http://milieuinterieur.fr"), align = "center"),
    
    tabsetPanel(type="tabs",
                tabPanel("Overview",
                         includeMarkdown("Introduction.md")),
                tabPanel("PCA View",
                         fluidPage(h2("Principal Component Analysis on the gene expression signature."),
                                   #p("This view ..."),
                                   h6("for more information, check the Instructions tab."),
                                   fluidRow(column(6, uiOutput("pcaFondationalStimuliSelectControl")),
                                            column(1, uiOutput("fig1AButton")),
                                            column(1, uiOutput("fig1BButton")),
                                            column(1, uiOutput("figS4AButton")),
                                            column(1, uiOutput("figS4BButton")),
                                            column(1, uiOutput("figS5Button"))),
                                   fluidRow(column(6,plotOutput("pcaPlots")), 
                                            column(6,plotOutput("pcaPlots2"))),
                                   fluidRow(column(5, uiOutput("pcaProjectedStimuliSelectControl"))),
                                   fluidRow(column(2, uiOutput("pcaGeneSetControl")),
                                            column(10,uiOutput("pcaGeneSetList"))))),
                tabPanel("Boxplot View",
                         h2("Gene expression levels across multiple stimuli."),
                         #p("This view ..."),
                         h6("for more information, check the Instructions tab."),
                         fluidRow(column(4, uiOutput("boxplotStimuliSelectControl")), 
                                  column(4, uiOutput("boxplotmRNAsSelectControl")), 
                                  column(1, uiOutput("fig2AButton")),
                                  column(1, uiOutput("fig2CButton")),
                                  column(1, uiOutput("fig3BButton"))), 
                         fluidRow(downloadButton('downloadmRNABoxplot', 'Download Plot'),
                                  column(2,checkboxInput("boxplot_muliple_labels", label = "Compact plots", value = FALSE))), 
                         plotOutput("mRNAsBoxplots")),
                tabPanel("Correlations",
                         h2("Gene expression correlations across stimuli."),
                         h6("for more information, check the Instructions tab."),
                         #p("This view ..."),
                         uiOutput("stimuliCorrelationmRNASelectControl"), 
                         plotOutput("stimuliCorrelationDendogram"), 
                         plotOutput("stimuliCorrelationMatrix")),
                tabPanel("Reference Values",
                         h2("Reference values table based on 25 healthy donors."),
                         h6("for more information, check the Instructions tab."),
                         fluidRow(column(4, selectInput("valueType", "Reference Value:",
                                                        c("Median expression" = "median", 
                                                          "Median fold change (versus Null)" = "fc",
                                                          "Coefficient of Variation (%)" = "cv", 
                                                          "q-value (paired t-test with Null)" = "qvalue")
                         )), column(5, uiOutput("refTableStimuliSelectControl")), column(3,numericInput("refTableMaxQvalue", "Highlight cells: define max paired t-test q-value ", 0.001, min = 0.0, max = 0.01, step = 0.001))),
                         DT::dataTableOutput("dataTable"), downloadButton('downloadReferenceValues', 'Download')),
                tabPanel("Instructions",includeMarkdown("Instructions.md")))
))