library(shiny)

shinyServer(function(input, output) {
  output$test <- renderText({
    input$test
  })
})

# vim: set ts=2 sw=2 foldmethod=marker expandtab:
