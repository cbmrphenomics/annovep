#!/usr/bin/env Rscript

requireNamespace("shiny")
requireNamespace("DT")

example_query <- paste(
  "# Select variants with a MODERATE/HIGH impact and less than 10% frequency in non-Finnish Europeans:",
  "(Func_impact = 'HIGH' OR Func_impact = 'MODERATE') AND gnomAD_NFE_AF < 0.1",
  sep = "\n"
)

ui <- pageWithSidebar(
  headerPanel("AnnoVEP"),
  sidebarPanel(
    # Authentication
    passwordInput("password", "Please enter password:"),
    uiOutput("db_errors"),
    hr(),
    conditionalPanel(
      condition = "output.show_ui && !output.database_errors",
      # View by Genes / Chromosome
      radioButtons("build", "Build:", c("Hg38" = "hg38", "Hg19 (liftover)" = "hg19"), inline = TRUE),
      radioButtons("select_by", "Select variants by:", c("Gene" = "genes", "Chromosome" = "chromosome"), inline = TRUE),
      conditionalPanel(
        condition = "input.select_by === 'genes'",
        selectizeInput("genes", NULL, choices = NULL, multiple = TRUE, options = list(maxOptions = 10, maxItems = 10))
      ),
      conditionalPanel(
        condition = "input.select_by === 'chromosome'",
        selectInput("chr", NULL, choices = NULL),
        numericInput("min_pos", "Start position (bp):", 1, min = 1, step = 1),
        numericInput("max_pos", "End position (bp; max 10k rows shown):", NULL, min = 1, step = 1)
      ),
      hr(),
      # Filters
      checkboxInput("require_pass", "Variants must PASS quality filtering", value = TRUE),
      selectInput("consequence",
        HTML("This consequence or worse (<a target='_blank' href='https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html'>definitions</a>)"),
        choices = c("Any consequence")
      ),
      numericInput("min_maf", "Minimum MAF (gnomAD):", 0, min = 0, max = 1, step = 0.00001),
      numericInput("max_maf", "Maximum MAF (gnomAD):", 1, min = 0, max = 1, step = 0.00001),
      textAreaInput("query", "Filters", value = "", placeholder = example_query, rows = 3),
      uiOutput("query_errors"),
      actionButton("btn_reset", "Reset filters"),
      hr(),
      # Visible columns
      selectizeInput("columns", "Visible columns", choices = c(), multiple = TRUE),
      hr(),
      # Export buttons
      downloadButton("btn_download", "Download")
    ),
    width = 3
  ),
  mainPanel(
    DT::dataTableOutput("table"),
    width = 9
  )
)
