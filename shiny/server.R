#!/usr/bin/env Rscript

library(DBI)
library(shiny)
library(DT)


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
    validate(need(expr, sprintf("Error in query near '%s': %s", token, message)))
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

  if (is.null(query)) {
    # Initial state
  } else if (open_brackets > 0) {
    fail(tail(tokens, n = 1), "unbalanced brackets")
  } else if (!is_state("rbracket", "value")) {
    fail(tail(tokens, n = 1), token, "partial query")
  }

  return(list(string = paste0(query, collapse = ""), params = params))
}


# Read-only connection to SQLite3 database
database <- Sys.getenv("ANNOVEP_DATABASE", "sqlite3.db")
conn <- dbConnect(RSQLite::SQLite(), database, flags = RSQLite::SQLITE_RO)

# FIXME: Better solution needed
password <- Sys.getenv("ANNOVEP_PASSWORD", ids::uuid())
cat("Password is", password, end = "\n")

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


server <- function(input, output, session) {
  userQuery <- reactiveValues(
  # The last valid user query
    valid_query = list(string = "", params = list()),
  # Error messages from parsing query, if any
    errors = NULL
  )

  observeEvent({
    input$query
    input$password
  }, {
    if (input$password == password) {
      tryCatch({
        query <- parse_query(input$query, symbols = columns)
        # Avoid spurious updates for white space changes
        if (!identical(query, userQuery$valid_query)) {
          userQuery$valid_query <- query
        }

        userQuery$errors <- NULL
      },
      error = function(cond) {
        userQuery$errors <- cond
      })
    } else {
      userQuery$valid_query <- list(string = "", params = list())
      userQuery$errors <- NULL
    }
  })

  # Returns the user-selected region of interest
  select_region <- function(input) {
    validate(need(input$password == password, "Password required"))

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

    user_query <- userQuery$valid_query
    if (nchar(user_query$string) > 0) {
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

  output$uiQueryErrors <- renderUI({
    if (!is.null(userQuery$errors)) {
      span(
        HTML("<h5 style='color: #AB0000; text-align: center;'>"),
        userQuery$errors,
        HTML("</h5>")
      )
    }
  })

  # Fill in dynamic list of consequence terms
  updateSelectInput(session, "consequence", choices = c("Any consequence", consequences))

  # Enable automatic reconnections
  session$allowReconnect("force")
}
