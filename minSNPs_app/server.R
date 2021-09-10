library(shiny)
library(shinyjs)
library(minSNPs)

shinyServer(function(input, output, session) {

  analysis_mode <- reactive({
    x <- input$metric
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(paste("Analysis Mode:", x))
  })

  process_mode <- reactiveValues()

  observeEvent({
      input$file1
    }, {
    inFile <- input$file1
    if (!is.null(inFile)) {
      analysis_sequences <- read_fasta(inFile$datapath)
      processed_sequences <- process_allele(analysis_sequences)
      process_mode[["processed_sequences"]] <- processed_sequences
      enable("run")
    } else {
      disable("run")
    }
  })

  to_listen <- reactive({
    list(process_mode[["processed_sequences"]], input$metric)
  })

  observeEvent(
      to_listen(), {
    metric <- input$metric
    sequences <- process_mode[["processed_sequences"]]
     if (metric == "percent" &&
        !is.null(sequences)) {
      updateSelectInput(session, "goi",
        choices = names(sequences$seqc)
      )
    } else {
      updateSelectInput(session, "goi",
        choices = list()
      )
    }
  })

  run_analysis <- eventReactive(input$run, {
    if (! is.null(process_mode[["processed_sequences"]])) {
      disable("run")
      shinyjs::show("status")
      process_mode[["results"]] <- find_optimised_snps(
        seqc = process_mode[["processed_sequences"]],
        metric = input$metric, goi = input$goi,
        accept_multiallelic = input$accept_multiallelic,
        included_positions = as.numeric(input$included_positions),
        excluded_positions = as.numeric(input$excluded_positions),
        number_of_result = input$number_of_result,
        max_depth = input$max_depth, bp = BiocParallel::MulticoreParam()
      )
      enable("run")
      shinyjs::hide("status")
      enable("download_result")
    }
  })

  output$status <- renderText({
    "runnning"
  })

  output$mode <- renderText({
    analysis_mode()
  })

  output$processed <- renderPrint({
    processed <- process_mode[["processed_sequences"]]
    if (!is.null(processed)) {
      print("Ignored positions:")
      print(paste(processed[["ignored_position"]], collapse = ","))
      print("Ignored isolates:")
      print(paste(processed[["ignored_allele"]], collapse = ","))
    }
  })

  output$download_result <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", as.character(Sys.time())),
        ".tsv", sep = "")
    },
    content = function(file) {
      output_result(process_mode[["results"]],
        view = "csv", file_name = file,
        seqc = process_mode[["processed_sequences"]])
    },
    contentType = "text/csv"
  )

  output$summary <- renderPrint({
    run_analysis()
    if (!is.null(process_mode[["results"]]) &&
      !is.null(process_mode[["processed_sequences"]])) {
        enable("download_result")
        output_result(process_mode[["results"]], view = "",
          seqc = process_mode[["processed_sequences"]])
    } else{
      "No Result"
    }
  })
})
