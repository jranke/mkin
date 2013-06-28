library(shiny)
shinyUI(pageWithSidebar(

  headerPanel("shiny GUI for mkin"),

  sidebarPanel(
    textInput("test", "Test:", "test input")
  ),

  mainPanel(
    textOutput("test")
  )
))
# vim: set ts=2 sw=2 foldmethod=marker expandtab:
