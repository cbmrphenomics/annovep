#!/usr/bin/env Rscript

# Enable backtrace on errors
options(error = function(...) {
  traceback(2)
  quit(status = 1)
})

has.values <- function(...) {
  for (value in list(...)) {
    if (is.null(value) || length(value) == 0) {
      return(FALSE)
    }
  }

  return(TRUE)
}


main <- function(args) {
  if (length(args) < 1) {
    cat("Usage:\n")
    cat("Rscript <sqlite.db>\n")
    return(1)
  }

  # TODO: More stringent requirements for length, etc.
  password <- Sys.getenv("ANNOVEP_PASSWORD", NA)
  if (is.na(password)) {
    cat("Please set environmental variable `ANNOVEP_PASSWORD` with desired password")
    return(1)
  }

  library(DBI)
  library(shiny)
  library(DT)

  # Read-only connection to SQLite3 database
  conn <- dbConnect(RSQLite::SQLite(), args[1], flags = RSQLite::SQLITE_RO)

  dbQuery <- function(string, ...) {
    return(dbGetQuery(conn, string, ...))
  }

  dbQueryVec <- function(string, ...) {
    return(unlist(dbQuery(string, ...), use.names = FALSE))
  }

  chroms <- dbQueryVec("SELECT DISTINCT [Chr] FROM [Annotations];")
  columns <- dbQueryVec("SELECT [Name] FROM [Columns];")
  consequences <- dbQueryVec("SELECT [Name] FROM [Consequences];")
  genes <- sort(dbQueryVec("SELECT [Name] FROM [Genes];"))

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

      radioButtons("filter", "VCF filters:", selected = "pass", c("PASS" = "pass", "Any value" = "any"), inline = TRUE),
      selectInput("consequence", "This consequence or worse", choices = c("Any consequence", consequences)),
      uiOutput("uiGene"),

      width = 3
    ),

    mainPanel(
      DT::dataTableOutput("table"),

      width = 9
    )
  )

  server <- function(input, output, session) {
    # Returns the user-selected region of interest
    select_region <- function(input) {
      if (input$select_by == "gene") {
        result <- dbQuery("SELECT * FROM [Genes] WHERE Name = :name ", params = list(name = input$gene))

        list(chr = result$Chr, minPos = result$Start, maxPos = result$End)
      } else {
        list(chr = input$chr, minPos = input$minPos, maxPos = maxPos())
      }
    }

    apply_filters <- function(input) {
      filters <- NULL
      if (input$filter == "pass") {
        filters <- "  AND Filters = 'PASS'"
      }

      consequence_idx <- match(input$consequence, consequences)
      if (!is.na(consequence_idx)) {
        filters <- c(filters, sprintf("  AND Func_most_significant <= %i", consequence_idx))
      }

      return(paste(c(filters, ""), collapse = "\n"))
    }

    apply_max_rows <- function(input) {
      if (input$select_by == "region") {
        return(sprintf("LIMIT %i", max(1, min(10000, input$maxRows))))
      } else {
        return("")
      }
    }

    create_sort_order <- function(table, column) {
      column_order <- sprintf("%s_order", column)
      order <- table[, column]

      table[, column_order] <- order
      table[, column] <- consequences[order]
      table[is.na(order), column_order] <- length(consequences) + 1

      return(table)
    }

    observeEvent({
      input$password
      input$select_by
    }, {
      visible_genes <- if (input$password == password) { genes } else { NULL }
      updateSelectizeInput(session, "gene", choices = visible_genes, server = TRUE)
    })


    min_max_position <- function(input, what) {
      if (has.values(input$chr)) {
        query <- sprintf("SELECT %s(Pos) FROM [Annotations] WHERE Chr = :chr%s", what, apply_filters(input))

        dbQueryVec(query, params = list(chr = input$chr))
      } else {
        1
      }
    }

    minPos <- reactive({ min_max_position(input, "MIN") })
    maxPos <- reactive({ min_max_position(input, "MAX") })

    output$uiChrom <- renderUI({
      visible_chroms <- if (input$password == password) { chroms } else { "" }
      selectInput("chr", NULL, choices = visible_chroms, selected = visible_chroms[1])
    })

    output$uiOffset <- renderUI({
      visible_min <- if (input$password == password) { minPos() } else { 1 }
      visible_max <- if (input$password == password) { maxPos() } else { 1 }
      numericInput("minPos", "Start position (bp):", visible_min, min = visible_min, max = visible_max, step = 1)
    })

    output$table <- DT::renderDataTable({
      validate(need(input$password == password, "Password required"))

      region <- select_region(input)
      if (has.values(region$chr, region$minPos, region$maxPos)) {
        query <- sprintf("
            SELECT *
            FROM   [Annotations]
            WHERE  Chr = :chr
              AND  Pos >= :min
              AND  Pos <= :max
            %s
            %s;",
            apply_filters(input),
            apply_max_rows(input))

        result <- dbQuery(query, params = list(chr = region$chr, min = region$minPos, max = region$maxPos))

        result$pid <- NULL # Not relevant
        if (nrow(result) > 0) {
          result <- create_sort_order(result, "Func_most_significant")
          result <- create_sort_order(result, "Func_least_significant")
          result <- create_sort_order(result, "Func_most_significant_canonical")
        }

        result
      }
    },
      selection = "single",
      options = list(
        scrollX = TRUE,
        scrollY = "80vh",
        pageLength = 100,
        lengthMenu = list(c(50, 100, 200, -1), c(50, 100, 200, "All")),
        columnDefs = list(
    # Order consequences by their rank using values from hidden columns
    # Note that while targets are 0-indexed, this includes an automatic column containing row numbers
          list(targets = match("Func_most_significant", columns), orderData = length(columns) + 1),
          list(targets = match("Func_least_significant", columns), orderData = length(columns) + 2),
          list(targets = match("Func_most_significant_canonical", columns), orderData = length(columns) + 3),
          list(targets = c(-1, -2, -3), visible = FALSE),
          list(targets = match("Info", columns), visible = FALSE)
        )
      )
    )


    # Enable automatic reconnections
    session$allowReconnect("force")
  }

  cat("starting server ...\n")

  runApp(shinyApp(ui, server, options = list(port = 56789)))

  return(0)
}


quit(status = main(commandArgs(TRUE)))
