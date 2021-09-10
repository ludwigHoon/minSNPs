library(shiny)
library(shinyjs)

shinyUI(
  pageWithSidebar(
    headerPanel("MinSNPs Analysis"),

    sidebarPanel(
      shinyjs::useShinyjs(),
      fileInput("file1", "FASTA FILE",
                accept = c(".fasta", ".fas", ".fa", ".txt")),
      selectInput("metric", "Analysis Mode:",
                  list("%-mode" = "percent",
                      "d-mode" = "simpson")),
      numericInput("number_of_result", "Number of result:", 5),
      numericInput("max_depth", "Max Depth:", 5),
      selectInput("goi", "Group of Interest:", choices = list(),
                  multiple = TRUE),
      selectizeInput("included_positions", "Included positions",
        choices = NULL, multiple = TRUE, options = list(create = TRUE)),
      selectizeInput("excluded_positions", "Excluded positions",
        choices = NULL, multiple = TRUE, options = list(create = TRUE)),
      checkboxInput("accept_multiallelic", "Accept Multiallelic SNP in GOI",
        FALSE),
      disabled(actionButton("run", "Start analysis")),
      disabled(downloadButton("download_result", "Download Result"))
    ),

    mainPanel(
      h3(textOutput("mode")),
      hidden(h4(textOutput("status"))),
      verbatimTextOutput("processed"),
      verbatimTextOutput("summary"),
    )
  )
)
