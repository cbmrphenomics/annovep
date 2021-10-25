#!/usr/bin/env Rscript

library(DBI)
library(shiny)
library(DT)


hasValues <- function(...) {
  for (value in list(...)) {
    if (is.null(value) || length(value) == 0) {
      return(FALSE)
    } else if (is.character(value) && max(nchar(value)) == 0) {
      return(FALSE)
    }
  }

  return(TRUE)
}

withDefault <- function(value, default) {
  if (hasValues(value)) {
    return(value)
  } else {
    return(default)
  }
}


parse_query <- function(value, symbols, special_values = list()) {
  # Parsing rules for simplified SQL "WHERE" conditions
  # In addition, the following terms are used
  #   logical: symbols AND or OR
  #   value:   number | string | symbol
  rules <- c(
    number = paste("\\b", flexo::re$number, "\\b", sep = ""),
    string = "\"[^\"]*\"|'[^']*'",
    symbol = "\\b\\w+\\b",
    operator = "==|=|!=|<>|<=|<|>=|>",
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
      if (name == "symbol" && toupper(token) == "LIKE") {
        name <- "operator"
        token <- "LIKE"
      }

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

        # Some values may be strings mapped onto numbers, etc.
        value <- withDefault(special_values[[tolower(token)]], token)

        key <- sprintf("up%i", length(params) + 1)
        query <- c(query, ":", key)
        params[key] <- value
      } else {
        fail(token, "expected value after operator")
      }
    } else if (is_state("value", "rbracket")) {
      if (name == "symbol") {
        symbol <- toupper(token)
        if (symbol == "AND" || symbol == "OR") {
          name <- "logical"
          query <- c(query, " ", symbol, " ")
        } else {
          fail(token, "unexpected symbol; expected AND or OR")
        }
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
    fail(tail(tokens, n = 1), "partial query")
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

chroms <- dbQueryVec("SELECT DISTINCT [Chr] FROM [Annotations] ORDER BY [Chr];")
columns <- dbQueryVec("SELECT [Name] FROM [Columns] ORDER BY [pid];")
consequences <- dbQuery("SELECT [pid], [Name] FROM [Consequences] ORDER BY [pid];")
genes <- dbQueryVec("SELECT [Name] FROM [Genes] ORDER BY [Name];")

default_columns <- c(
  "Chr",
  "Pos",
  "Ref",
  "Alt",
  "Freq",
  "Ancestral_allele",
  "Func_most_significant",
  "Func_most_significant_canonical",
  "Fund_gene_symbol",
  "Func_protein_position",
  "Func_amino_acids",
  "dbSNP_ids"
)

# Consequences are foreign keys/ranks to allow ordering comparisons
special_values <- list()
special_values[consequences$Name] <- consequences$pid


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
        query <- parse_query(input$query, symbols = columns, special_values = special_values)
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
    consequence_pid <- consequences$pid[consequences$Name == input$consequence]
    if (hasValues(consequence_pid)) {
      query <- c(query, sprintf("  AND Func_most_significant >= %i", consequence_pid))
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

  create_sort_order <- function(table, visible_columns, column) {
    if (column %in% visible_columns) {
      column_order <- sprintf("%s_order", column)
      order <- table[, column]

      table[, column_order] <- -order
      table[, column] <- consequences$Name[match(order, consequences$pid)]
      table[is.na(order), column_order] <- 1
    }

    return(table)
  }

  observeEvent({
    input$password
  }, {
    if (input$password == password) {
      visible_genes <- genes
      selected_gene <- withDefault(input$gene, genes[1])
      visible_chroms <- chroms
      selected_chrom <- withDefault(input$chr, chroms[1])
      visible_columns <- columns
      selected_columns <- withDefault(input$columns, default_columns)
    } else {
      visible_genes <- NULL
      selected_gene <- NULL
      visible_chroms <- NULL
      selected_chrom <- NULL
      visible_columns <- NULL
      selected_columns <- NULL
    }

    updateSelectInput(session, "chr", selected = selected_chrom, choices = visible_chroms)
    updateSelectizeInput(session, "gene", selected = selected_gene, choices = visible_genes, server = TRUE)
    updateSelectizeInput(session, "columns", selected = selected_columns, choices = visible_columns)
  })

  maxPos <- reactive({
    if (hasValues(input$chr)) {
      query <- sprintf("SELECT MAX(Pos) FROM [Annotations] WHERE Chr = :chr")
      query <- build_query(input, query, params = list(chr = input$chr))

      dbQueryVec(query$string, params = query$params)
    } else {
      1
    }
  })

  observeEvent({
    input$password
    maxPos
  }, {
    max_value <- if (input$password == password) { maxPos() } else { 1 }
    value <- ifelse(hasValues(input$minPos), max(1, min(max_value, input$minPos)), 1)

    updateNumericInput(session, "minPos", value = value, max = max_value)
  })

  data <- reactive({
    validate(need(input$password == password, "Password required"))

    region <- select_region(input)
    # Ensure that only valid column names are used
    visible_columns <- subset(columns, columns %in% input$columns)
    if (hasValues(region$chr, region$minPos, region$maxPos, visible_columns)) {
      params <- list(chr = region$chr, min = region$minPos, max = region$maxPos)
      query <- c(
          sprintf("SELECT %s", paste(sprintf("[%s]", visible_columns), collapse=", ")),
          "FROM   [Annotations]",
          "WHERE  Chr = :chr",
          "  AND  Pos >= :min",
          "  AND  Pos <= :max"
      )

      query <- build_query(input, query, params)
      result <- dbQuery(query$string, params = query$params)

      result <- create_sort_order(result, visible_columns, "Func_most_significant")
      result <- create_sort_order(result, visible_columns, "Func_least_significant")
      result <- create_sort_order(result, visible_columns, "Func_most_significant_canonical")

      result
    }
  })

  observeEvent({
    data
    input$columns
  }, {
    consequence_columns <- c(
      "Func_most_significant",
      "Func_least_significant",
      "Func_most_significant_canonical"
    )

    # Ensure that only valid column names are used
    visible_columns <- subset(columns, columns %in% input$columns)
    visible_consequences <- subset(consequence_columns, consequence_columns %in% visible_columns)

    coldefs <- list()
    offset <- length(visible_columns)
    hidden_columns <- numeric(0)
    for (name in visible_consequences) {
      coldefs <- c(coldefs, list(list(targets = match(name, visible_columns), orderData = offset + 1)))
        hidden_columns <- c(hidden_columns, offset + 1)
        offset <- offset + 1
    }

    if (length(hidden_columns) > 0) {
    coldefs <- c(coldefs, list(list(targets = hidden_columns, visible = FALSE)))
    }

    output$table <- DT::renderDataTable(
      data(),
      selection = "single",
      server = TRUE,
      options = list(
        scrollX = TRUE,
        scrollY = "80vh",
        pageLength = 100,
        lengthMenu = list(c(50, 100, 200, -1), c(50, 100, 200, "All")),
        columnDefs = coldefs
      )
    )
  })

  output$uiQueryErrors <- renderUI({
    if (hasValues(userQuery$errors)) {
      span(
        HTML("<h5 style='color: #AB0000; text-align: center;'>"),
        userQuery$errors,
        HTML("</h5>")
      )
    }
  })

  # Fill in dynamic list of consequence terms
  updateSelectInput(session, "consequence", choices = c("Any consequence", rev(consequences$Name)))

  # Enable automatic reconnections
  session$allowReconnect("force")
}
