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
      uiOutput("uiChrom"),
      uiOutput("uiOffset"),
      numericInput("maxRows", "Max rows:", 10000, min = 1, max = 10000, step = 1)
    ),

    hr(),

    textAreaInput("query", "Filters", value = "Filters = PASS", rows = 3),
    uiOutput("uiQueryErrors"),
    selectInput("consequence", "This consequence or worse", choices = c("Any consequence")),

    width = 3
  ),

  mainPanel(
    DT::dataTableOutput("table"),

    width = 9
  )
)
