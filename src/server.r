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


parse_query <- function(value, symbols) {
  # Parsing rules for simplified SQL "WHERE" conditions
  rules <- c(
    number = paste("\\b", flexo::re$number, "\\b", sep = ""),
    logical = "\\b[aA][nN][dD]\\b|\\b[oO][rR]\\b",
    string = "\"[^\"]*\"|'[^']*'",
    symbol = "\\b\\w+\\b",
    operator = "==|=|!=|<>|<|<=|>|>=",
    lbracket = "\\(",
    rbracket = "\\)",
    whitespace = "\\s+"
  )

  symbols_lc <- tolower(symbols)
  tokens <- flexo::lex(value, rules)
  names <- names(tokens)
  query <- NULL
  params <- list()
  open_brackets <- 0

  # Safe starting state
  state <- c("whitespace", "logical")

  is_state <- function(...) {
    for (name in list(...)) {
      if (state[1] == name || (state[1] == "whitespace" && state[2] == name)) {
        return(TRUE)
      }
    }

    return(FALSE)
  }

  check <- function(expr, token, message) {
    validate(need(expr, sprintf("Error in query near '%s':\n  %s", token, message)))
  }

  fail <- function(token, message) { check(FALSE, token, message) }

  for (i in seq_along(tokens)) {
    token <- tokens[i]
    name <- names[i]

    if (name == "whitespace") {
      # always allowed
    } else if (is_state("logical", "lbracket")) {
      if (name == "lbracket") {
        open_brackets <- open_brackets + 1
        query <- c(query, "(")
      } else if (name == "string" || name == "symbol") {
        if (name == "string") {
          # Remove quotes and promote to symbol
          token <- substr(token, 2, nchar(token) - 1)
          name <- "symbol"
        }

        index <- match(tolower(token), symbols_lc)
        check(!is.na(index), token, "not a valid column name")

        query <- c(query, symbols[index])
      } else {
        fail(token, "expected brackets or a column name")
      }
    } else if (is_state("symbol")) {
      check(name == "operator", token, "expected operator")
      query <- c(query, " ", token, " ")
    } else if (is_state("operator")) {
      if (name == "string" || name == "symbol" || name == "number") {
        if (name == "string") {
          token <- substr(token, 2, nchar(token) - 1)
        } else if (name == "number") {
          token <- as.numeric(token)
        }

        # Downgrade everything to "value" to simplify logic
        name <- "value"

        key <- sprintf("up%i", length(params) + 1)
        query <- c(query, ":", key)
        params[key] <- token
      } else {
        fail(token, "expected value after operator")
      }
    } else if (is_state("value", "rbracket")) {
      if (name == "logical") {
        query <- c(query, " ", toupper(token), " ")
      } else if (name == "rbracket") {
        check(open_brackets >= 1, token, "unbalanced brackets")
        open_brackets <- open_brackets - 1
        query <- c(query, ")")
      } else {
        fail(token, "expected AND/OR or bracket")
      }
    } else {
      fail(token, "malformed query")
    }

    state <- c(name, state[1])
  }

  if (open_brackets > 0) {
    fail(tail(tokens, n = 1), "unbalanced brackets")
  } else if (!is_state("rbracket", "value")) {
    fail(token, "partial query")
  }

  return(list(input = value, string = paste0(query, collapse = ""), params = params))
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

      textAreaInput("query", "Filters", value = "Filters = PASS", rows = 3),
      selectInput("consequence", "This consequence or worse", choices = c("Any consequence", consequences)),

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

    build_query <- function(input, query, params) {
      consequence_idx <- match(input$consequence, consequences)
      if (!is.na(consequence_idx)) {
        query <- c(query, sprintf("  AND Func_most_significant <= %i", consequence_idx))
      }

      user_query <- parse_query(input$query, symbols = columns)
      if (!is.null(user_query$string)) {
        query <- c(query, paste("AND", user_query$string, sep = " "))
        params <- c(params, user_query$params)
      }

      if (input$select_by == "region") {
        query <- c(query, sprintf("LIMIT %i", max(1, min(10000, input$maxRows))))
      }

      return(list(string = paste(c(query, ";"), collapse = "\n"), params = params))
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

    maxPos <- reactive({
      if (has.values(input$chr)) {
        query <- sprintf("SELECT MAX(Pos) FROM [Annotations] WHERE Chr = :chr")
        query <- build_query(input, query, params = list(chr = input$chr))

        dbQueryVec(query$string, params = query$params)
      } else {
        1
      }
    })

    output$uiChrom <- renderUI({
      visible_chroms <- if (input$password == password) { chroms } else { "" }
      selectInput("chr", NULL, choices = visible_chroms, selected = visible_chroms[1])
    })

    output$uiOffset <- renderUI({
      visible_min <- isolate(input$minPos)
      visible_min <- ifelse(is.null(visible_min), 1, visible_min)
      visible_max <- if (input$password == password) { maxPos() } else { 1 }
      numericInput("minPos", "Start position (bp):", visible_min, min = 1, max = visible_max, step = 1)
    })

    output$table <- DT::renderDataTable({
      validate(need(input$password == password, "Password required"))

      region <- select_region(input)
      if (has.values(region$chr, region$minPos, region$maxPos)) {
        params <- list(chr = region$chr, min = region$minPos, max = region$maxPos)
        query <- c(
            "SELECT *",
            "FROM   [Annotations]",
            "WHERE  Chr = :chr",
            "  AND  Pos >= :min",
            "  AND  Pos <= :max"
        )

        query <- build_query(input, query, params)
        result <- dbQuery(query$string, params = query$params)

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
      server = TRUE,
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
