library(shiny)

ui <- fluidPage(
  titlePanel("Allele Specific Expression Data"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("geneList", 
                 "Choose a subset of genes",
                 choices = c("ASE" = "aseGenesBeta",
                             "DE" = "DEGenes", 
                             "ASE and DE" = "ASEDEGenes",
                             "All Genes" = "genes")),
    
      uiOutput("geneControls")
    ),
                 
    mainPanel(
      textOutput("selectedGene"),
      plotOutput("barPlot"),
      plotOutput("sleuthPlot"),
      textOutput("GOinfo")
      )
    )
  )
