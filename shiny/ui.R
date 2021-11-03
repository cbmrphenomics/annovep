#!/usr/bin/env Rscript

library(shiny)
library(DT)

ui <- pageWithSidebar(
  headerPanel("AnnoVEP"),

  sidebarPanel(
    passwordInput("password", "Password"),

    hr(),

    radioButtons("select_by", "Select by:", c("Gene" = "gene", "Region" = "region"), inline = TRUE),

    conditionalPanel(
      condition = "input.select_by === 'gene'",
      selectizeInput("gene", NULL, choices = NULL, options = list(maxOptions = 10))
    ),

    conditionalPanel(
      condition = "input.select_by === 'region'",
      selectInput("chr", NULL, choices = NULL),
      numericInput("minPos", "Start position (bp):", 1, min = 1, step = 1),
      numericInput("maxPos", "End position (bp; max 10k rows shown):", NULL, min = 1, step = 1)
    ),

    hr(),

    textAreaInput("query", "Filters", value = "Filters = PASS", rows = 3),
    uiOutput("uiQueryErrors"),
    selectInput("consequence", "This consequence or worse", choices = c("Any consequence")),

    hr(),

    selectizeInput("columns", "Visible columns", choices = c(), multiple = TRUE),

    width = 3
  ),

  mainPanel(
    DT::dataTableOutput("table"),

    width = 9
  )
)
