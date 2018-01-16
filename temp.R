  
  library(shiny)
  
  
  ui = bootstrapPage(
      numericInput('n', 'Number of bins', 100),
      plotOutput('plot')
    )
  server = function(input, output) {
      output$plot <- renderPlot({ hist(dat, breaks = input$n) })
    }

