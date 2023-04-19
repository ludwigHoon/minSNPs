library(shiny)
source("private_dependency.R")

shiny::runApp(".", port = 3838, launch.browser = TRUE)
