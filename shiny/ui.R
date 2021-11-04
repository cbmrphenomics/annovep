#!/usr/bin/env Rscript

requireNamespace("shiny")
requireNamespace("DT")

example_query <- paste(
  "# Select variants with a MODERATE/HIGH impact and less than 10% frequency in non-Finnish Europeans:",
  "(Func_impact = 'HIGH' OR Func_impact = 'MODERATE') AND gnomAD_NFE_AF < 0.1",
  sep="\n"
)

ui <- pageWithSidebar(
  headerPanel("AnnoVEP"),
  sidebarPanel(
    # Authentication
    passwordInput("password", "Password"),
    hr(),
    # View by Gene / region
    radioButtons("select_by", "Select by:", c("Gene" = "gene", "Region" = "region"), inline = TRUE),
    conditionalPanel(
      condition = "input.select_by === 'gene'",
      selectizeInput("gene", NULL, choices = NULL, options = list(maxOptions = 10))
    ),
    conditionalPanel(
      condition = "input.select_by === 'region'",
      selectInput("chr", NULL, choices = NULL),
      numericInput("min_pos", "Start position (bp):", 1, min = 1, step = 1),
      numericInput("max_pos", "End position (bp; max 10k rows shown):", NULL, min = 1, step = 1)
    ),
    hr(),
    # Filters
    checkboxInput("require_pass", "Variants must PASS quality filtering", value = TRUE),
    selectInput("consequence", "This consequence or worse", choices = c("Any consequence")),
    numericInput("min_maf", "Minimum MAF (gnomAD):", 0, min = 0, max = 1, step = 0.00001),
    numericInput("max_maf", "Maximum MAF (gnomAD):", 1, min = 0, max = 1, step = 0.00001),
    textAreaInput("query", "Filters", value = "", placeholder = example_query, rows = 3),
    uiOutput("query_errors"),
    actionButton("btn_reset", "Reset filters"),
    hr(),
    # Visible columns
    selectizeInput("columns", "Visible columns", choices = c(), multiple = TRUE),
    width = 3
  ),
  mainPanel(
    DT::dataTableOutput("table"),
    width = 9
  )
)
