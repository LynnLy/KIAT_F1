library(shiny)
library(ggplot2)
library(reshape2)
library(sleuth)
library(GO.db)

load("data/BngoList.Rdata")

# plot_bootstrap is taken from the sleuth shiny app
plot_bootstrap <- function(obj,
                           target_id,
                           units = 'est_counts',
                           color_by = setdiff(colnames(obj$sample_to_covariates), 'sample'),
                           x_axis_angle = 50,
                           divide_groups = TRUE
) {
  #units <- check_quant_mode(obj, units)
  
  df <- get_bootstrap_summary(obj, target_id, units)
  
  p <- ggplot(df, aes(x = sample, ymin = min, lower = lower, middle = mid,
                      upper = upper, ymax = max))
  p <- p + geom_boxplot(stat = "identity", aes_string(fill = color_by))
  p <- p + theme(axis.text.x = element_text(angle = x_axis_angle, hjust = 1))
  p <- p + ylab(units)
  if (divide_groups) {
    p <- p + facet_wrap(color_by, switch = 'x', scales = 'free_x')
  }
  
  p
}

shinyServer(function(input, output) {
  
  output$barPlot <- renderPlot({ 
    
    ggplot(visualDataTidy[visualDataTidy$GeneID == input$gene,]) + 
      geom_bar(aes(x = allele, y = reads, fill = allele), stat = "identity") + facet_wrap(~cultivar) +
      ylab("reads") 
    
  })
  
  output$geneControls <- renderUI({
    genesToUse <- get(input$geneList)
    selectizeInput("gene",
                   "Choose a gene",
                   choices=genesToUse)
  })
  
  output$sleuthPlot <- renderPlot({
    plot_bootstrap(so, input$gene, color_by = "cultivar")
  })
 
  output$selectedGene <- renderText({
    paste(input$gene)
  })
  
  output$GOinfo <- renderPrint({
    goTerms <- BngoList[[input$gene]]
    info <- vector(mode = "character", length = length(goTerms))
    info2 <- vector(mode = "character", length = length(goTerms))
    for(i in 1:length(goTerms)) {
      info[i] <- Term(GOTERM[[goTerms[i]]])
      info2[i] <- Definition(GOTERM[[goTerms[i]]])
    }
    info <- cbind(info, info2)
    info
    
  })
  
})
