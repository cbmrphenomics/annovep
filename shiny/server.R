#!/usr/bin/env Rscript

defaults <- list(
  password = paste0(sample(c(letters, LETTERS, 0:9), 16), collapse = ""),
  database = "sqlite3.db",
  build = "hg38",
  chrom = NULL,
  genes = NULL,
  columns = c(
    "Chr",
    "Pos",
    "Ref",
    "Alt",
    "Freq",
    "Ancestral_allele",
    "Genes_overlapping",
    "Func_most_significant",
    "Func_most_significant_canonical",
    "Fund_gene_symbol",
    "Func_protein_position",
    "Func_amino_acids",
    "dbSNP_ids",
    "gnomAD_min",
    "gnomAD_max"
  )
)

########################################################################################
## Logging functions

if (isatty(stdout())) {
  CLI_GREEN <- "\x1b[32m"
  CLI_GREY <- "\x1b[1;30m"
  CLI_RED <- "\x1b[31m"
  CLI_END <- "\x1b[0m"
} else {
  CLI_GREEN <- ""
  CLI_GREY <- ""
  CLI_RED <- ""
  CLI_END <- ""
}

log <- function(color, level, ...) {
  time <- strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(..., sep = "")

  cat(CLI_GREEN, time, CLI_END, " ", color, level, CLI_END, " ", msg, "\n", sep = "")

  return(invisible(msg))
}

info <- function(...) {
  return(log(CLI_GREY, "INFO", ...))
}

error <- function(...) {
  return(log(CLI_RED, "ERROR", ...))
}


########################################################################################
## Debugging functions


require_value <- function(check, min_len = 1, max_len = 1) {
  assert <- function(context, ..., optional = FALSE) {
    expr <- deparse(bquote(...))

    if (!is.character(context)) {
      stop(error("[?] context of assert is not a string: ", context))
    } else if (length(list(...)) == 0) {
      stop(error("[", context, "] empty expression"))
    }

    if (is.null(...)) {
      if (!optional) {
        stop(error("[", context, "] ", expr, " is NULL"))
      }
    } else if (!check(...)) {
      stop(error("[", context, "] ", expr, " has wrong type: ", typeof(...)))
    } else if (length(...) < min_len || length(...) > max_len) {
      stop(error("[", context, "] ", expr, " has wrong length (", length(...), ")"))
    }
  }

  return(assert)
}

require_bool <- require_value(is.logical)
require_num <- require_value(is.numeric)
require_str <- require_value(is.character)
require_nums <- require_value(is.numeric, min_len = 0, max_len = Inf)
require_strs <- require_value(is.character, min_len = 0, max_len = Inf)
require_list <- require_value(is.list, min_len = 0, max_len = Inf)

########################################################################################
## Initial setup

# Source optional R-file with user settings and custom library calls
load_settings <- function(settings) {
  if (file.exists("settings.R")) {
    info("Loading settings from 'settings.R'")
    source("settings.R", local = TRUE)
  } else {
    info("No 'settings.R' file found; using default settings")
  }

  require_list("user settings", settings)
  require_str("user settings", settings$password)
  require_str("user settings", settings$database)
  require_str("user settings", settings$genes, optional = TRUE)
  require_str("user settings", settings$build, optional = TRUE)
  require_str("user settings", settings$chrom, optional = TRUE)
  require_strs("user settings", settings$columns, optional = TRUE)

  return(settings)
}


settings <- load_settings(defaults)

info("Checking for required packages")
requireNamespace("DBI")
requireNamespace("DT")
requireNamespace("flexo")
requireNamespace("RSQLite")
requireNamespace("shiny")


has_values <- function(...) {
  for (value in list(...)) {
    if (is.null(value) || length(value) == 0) {
      return(FALSE)
    } else if (is.character(value) && max(nchar(value)) == 0) {
      return(FALSE)
    }
  }

  return(TRUE)
}

with_default <- function(value, default) {
  if (has_values(value)) {
    return(value)
  } else {
    return(default)
  }
}


parse_query <- function(value, symbols, special_symbols = list(), special_values = list()) {
  require_str("parse_query", value)
  require_strs("parse_query", symbols)
  require_list("parse_query", special_symbols)
  require_list("parse_query", special_values)

  # Parsing rules for simplified SQL "WHERE" conditions
  # In addition, the following terms are used
  #   logical: symbols AND or OR
  #   value:   number | string | symbol
  rules <- c(
    comment = "#.*",
    number = paste("\\b", flexo::re$number, "\\b", sep = ""),
    string = "\"[^\"]*\"|'[^']*'",
    symbol = "\\b\\w+\\b",
    operator = "==|=|!=|<>|<=|<|>=|>",
    lbracket = "\\(",
    rbracket = "\\)",
    whitespace = "\\s+"
  )
  symbols <- c("NULL", symbols)
  symbols_lc <- tolower(symbols)
  tokens <- flexo::lex(value, rules)
  names <- names(tokens)
  query <- NULL
  params <- list()
  open_brackets <- 0

  # Safe starting state
  state <- "logical"

  is_state <- function(...) {
    for (name in list(...)) {
      if (state == name) {
        return(TRUE)
      }
    }

    return(FALSE)
  }

  check <- function(expr, token, message) {
    shiny::validate(shiny::need(expr, sprintf("Error in query near '%s': %s", token, message)))
  }

  fail <- function(token, message) {
    check(FALSE, token, message)
  }

  for (i in seq_along(tokens)) {
    token <- tokens[i]
    name <- names[i]

    if (name == "whitespace" || name == "comment") {
      name <- state
    } else if (is_state("logical", "lbracket")) {
      if (name == "lbracket") {
        open_brackets <- open_brackets + 1
        query <- c(query, "(")
      } else if (name == "symbol" && toupper(token) == "NOT" && is_state("logical")) {
        # NOT is treated as a continuation of the current logical operatorq
        name <- "logical"

        query <- c(query, " NOT ")
      } else if (name == "string" || name == "symbol") {
        if (name == "string") {
          # Remove quotes and promote to symbol
          token <- substr(token, 2, nchar(token) - 1)
          name <- "symbol"
        }

        # Some values may be strings mapped onto numbers, etc.
        value <- with_default(special_symbols[[tolower(token)]], token)
        index <- match(tolower(value), symbols_lc)
        check(!is.na(index), token, "not a valid column name")

        query <- c(query, symbols[index])
      } else {
        fail(token, "expected brackets or a column name")
      }
    } else if (is_state("symbol")) {
      if (name == "symbol" && toupper(token) %in% c("LIKE", "IS")) {
        name <- "operator"
        token <- toupper(token)
      }

      check(name == "operator", token, "expected operator")
      query <- c(query, " ", token, " ")
    } else if (is_state("operator")) {
      if (name == "symbol") {
        index <- match(tolower(token), symbols_lc)
        check(!is.na(index), token, "not a valid name")

        query <- c(query, symbols[index])
      } else if (name == "string" || name == "number") {
        if (name == "string") {
          token <- substr(token, 2, nchar(token) - 1)
        } else if (name == "number") {
          token <- as.numeric(token)
        }

        # Some values may be strings mapped onto numbers, etc.
        value <- with_default(special_values[[tolower(token)]], token)

        key <- sprintf("up%i", length(params) + 1)
        query <- c(query, ":", key)
        params[key] <- value
      } else {
        fail(token, "expected value after operator")
      }

      # Downgrade everything to "value" to simplify logic
      name <- "value"
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

    state <- name
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


# Attempt to open the user-provided database
database <- tryCatch(
  {
    # Read-only connection to SQLite3 database
    conn <- DBI::dbConnect(RSQLite::SQLite(), settings$database, flags = RSQLite::SQLITE_RO)

    query <- function(string, ...) {
      return(DBI::dbGetQuery(conn, string, ...))
    }

    query_vec <- function(string, ...) {
      return(unlist(query(string, ...), use.names = FALSE))
    }

    list(
      conn = conn,
      errors = NULL,
      query = query,
      query_vec = query_vec,
      chroms_hg19 = query_vec("SELECT DISTINCT [Name] FROM [Contigs] WHERE [Build] = 'hg19' ORDER BY [pk];"),
      chroms_hg38 = query_vec("SELECT DISTINCT [Name] FROM [Contigs] WHERE [Build] = 'hg38' ORDER BY [pk];"),
      columns = query_vec("SELECT [Name] FROM [Columns] ORDER BY [pk];"),
      columns_info = query("SELECT [Name], [Description] FROM [Columns] ORDER BY [pk];"),
      consequences = query("SELECT [pk], [Name] FROM [Consequences] ORDER BY [pk];"),
      genes = query_vec("SELECT [Name] FROM [Genes] ORDER BY [Name];"),
      has_json = length(query_vec("SELECT 1 FROM [JSON];")) > 0
    )
  },
  error = function(e) {
    error("Error opening/querying database: ", e)

    list(
      conn = NULL,
      errors = as.character(e),
      query = function(...) (stop(e)),
      query_vec = function(...) (stop(e)),
      chroms_hg19 = character(),
      chroms_hg38 = character(),
      columns = character(),
      columns_info = data.frame(Name = character(), Description = character(), stringsAsFactors = FALSE),
      consequences = data.frame(pk = integer(), Name = character(), stringsAsFactors = FALSE),
      genes = character(),
      has_json = FALSE
    )
  }
)

require_strs("database", database$errors, optional = TRUE)
require_strs("database", database$chroms_hg19)
require_strs("database", database$chroms_hg38)
require_strs("database", database$columns)
require_strs("database", database$columns_info$Name)
require_strs("database", database$columns_info$Description)
require_nums("database", database$consequences$pk)
require_strs("database", database$consequences$Name)
require_strs("database", database$genes)
require_bool("database", database$has_json)

# Consequences are foreign keys/ranks to allow ordering comparisons
special_values <- list()
special_values[database$consequences$Name] <- database$consequences$pk

# Fill out default values not set by the user
if (is.null(settings$genes)) {
  settings$genes <- database$genes[1]
}

if (is.null(settings$build)) {
  settings$build <- "hg38"
}


if (is.null(settings$chrom)) {
  if (settings$build == "hg19") {
    settings$chrom <- database$chroms_hg19[1]
  } else {
    settings$chrom <- database$chroms_hg38[1]
  }
}

server <- function(input, output, session) {
  user_filter <- reactiveValues(
    # The last valid user query
    query = list(string = "", params = list()),
    # Error messages from parsing query, if any
    errors = NULL
  )

  # `input$build` is used in queries, so for safety the user value is not used directly
  hg_build <- function(input) {
    if (is.null(input$build) || input$build != "hg19") {
      return("hg38")
    } else {
      return("hg19")
    }
  }

  observeEvent(
    {
      input$query
      input$password
      input$build
    },
    {
      if (input$password == settings$password) {
        # Allow user queries using Chr/Pos instead of Hg38_pos/Hg38_chr
        if (hg_build(input) == "hg38") {
          special_symbols <- list("chr" = "hg38_chr", "pos" = "hg38_pos")
        } else {
          special_symbols <- list("chr" = "hg19_chr", "pos" = "hg19_pos")
        }

        tryCatch(
          {
            # Parse user-query
            query <- parse_query(
              input$query,
              symbols = database$columns,
              special_symbols = special_symbols,
              special_values = special_values
            )

            # Avoid spurious updates for white space changes
            if (!identical(query, user_filter$query)) {
              user_filter$query <- query
            }

            user_filter$errors <- NULL
          },
          error = function(cond) {
            user_filter$errors <- cond
          }
        )
      } else {
        user_filter$query <- list(string = "", params = list())
        user_filter$errors <- NULL
      }
    }
  )

  # Returns the user-selected regions of interest
  select_regions <- function(input) {
    shiny::validate(shiny::need(input$password == settings$password, "Password required"))

    if (input$select_by == "genes") {
      # The number of visible genes is hard-capped at 10 to prevent performance problems
      params <- list(genes = head(input$genes, n = 10))
      result <- database$query("SELECT * FROM [Genes] WHERE Name IN (:genes)", params = params)

      # For simplicity regions are always given in terms of hg38
      list(chr = result$Hg38_chr, min_pos = result$Hg38_start, max_pos = result$Hg38_end)
    } else {
      list(chr = input$chr, min_pos = input$min_pos, max_pos = input$max_pos)
    }
  }

  build_query <- function(input, region, params, columns) {
    build <- hg_build(input)
    col_chr <- sprintf("%s_chr", build)
    col_pos <- sprintf("%s_pos", build)

    ####################################################################################
    # SELECT
    query <- "SELECT "

    # Chr is a column created on demand based on the user-selected build
    if ("Chr" %in% columns) {
      query <- c(query, sprintf("[%s] AS [Chr], ", col_chr))
    }

    # Pos is a column created on demand based on the user-selected build
    if ("Pos" %in% columns) {
      query <- c(query, sprintf("[%s] AS [Pos], ", col_pos))
    }

    # Chr and Pos do not correspond to actual columns (see above) and have to be removed
    db_columns <- setdiff(columns, c("Chr", "Pos"))
    # The hg38 coordinates are required for some functionality, so include if missing
    db_columns <- c(db_columns, setdiff(c("Hg38_chr", "Hg38_pos"), db_columns))

    query <- c(query, paste(sprintf("[%s]", db_columns), collapse = ", "))

    ####################################################################################
    # FROM
    query <- c(query, "FROM [Annotations]")

    ####################################################################################
    # WHERE

    # Select rows in regions of interest (one for contigs, one or more for genes)
    params_chr <- as.list(region$chr)
    params_min <- as.list(region$min_pos)
    params_max <- as.list(region$max_pos)

    names(params_chr) <- sprintf("chr%i", seq_along(params_chr))
    names(params_min) <- sprintf("min%i", seq_along(params_min))
    names(params_max) <- sprintf("max%i", seq_along(params_max))

    if (is.numeric(region$max_pos)) {
      params <- c(params_chr, params_min, params_max)
      where <- sprintf("(%s = :%s AND %s >= :%s AND %s <= :%s)", col_chr, names(params_chr), col_pos, names(params_min), col_pos, names(params_max))
    } else {
      params <- c(params_chr, params_min)
      where <- sprintf("(%s = :%s AND %s >= :%s)", col_chr, names(params_chr), col_pos, names(params_min))
    }

    where <- sprintf("(%s)", paste(where, collapse = " OR "))

    # Exclude rows that do not have positions on the selected build
    where <- c(where, sprintf("AND %s IS NOT NULL", col_chr))

    # Map consequence term onto numerical ID
    consequence_pk <- database$consequences$pk[database$consequences$Name == input$consequence]
    if (has_values(consequence_pk)) {
      where <- c(where, sprintf("  AND Func_most_significant >= %i", consequence_pk))
    }

    # Whether or not the variant passed QC filters
    if (input$require_pass) {
      where <- c(where, "  AND Filters = 'PASS'")
    }

    # gnomAD min/max minor allele frequencies
    where <- c(where, "  AND gnomAD_min >= :gmin AND gnomAD_max <= :gmax")
    params <- c(params, list(gmin = input$min_maf, gmax = input$max_maf))

    # Optional user-provided filter
    user_query <- user_filter$query
    if (nchar(user_query$string) > 0) {
      where <- c(where, paste("AND", user_query$string, sep = " "))
      params <- c(params, user_query$params)
    }

    query <- c(query, "WHERE", where)

    ####################################################################################
    # ORDER BY
    query <- c(query, sprintf("ORDER BY %s, %s", col_chr, col_pos))

    ####################################################################################
    # LIMIT
    if (input$select_by == "contig") {
      query <- c(query, sprintf("LIMIT %i", max(1, min(10000, input$maxRows))))
    }

    return(list(
      string = paste(c(query, ";"), collapse = "\n"),
      params = params,
      columns = columns
    ))
  }

  create_sort_order <- function(table, visible_columns, column) {
    if (column %in% visible_columns) {
      column_order <- sprintf("%s_order", column)

      order <- table[, column]
      table[, column_order] <- ifelse(is.na(order), 1, order * -1)
      table[, column] <- database$consequences$Name[match(order, database$consequences$pk)]
    }

    return(table)
  }

  observeEvent(
    {
      input$password
    },
    {
      if (input$password == settings$password) {
        visible_genes <- database$genes
        selected_gene <- with_default(input$genes, settings$genes)
        visible_chroms <- database$chroms_hg38
        selected_chrom <- with_default(input$chr, settings$chrom)
        visible_columns <- c("Chr", "Pos", database$columns)
        selected_columns <- with_default(input$columns, settings$columns)

        updateTextInput(session, "password", "✔️ Password accepted!")
      } else {
        visible_genes <- NULL
        selected_gene <- NULL
        visible_chroms <- NULL
        selected_chrom <- NULL
        visible_columns <- NULL
        selected_columns <- NULL

        if (nchar(input$password) > 0) {
          updateTextInput(session, "password", "❌ Wrong password:")
        } else {
          updateTextInput(session, "password", "Please enter password:")
        }
      }

      shiny::updateRadioButtons(session, "build", selected = settings$build)
      shiny::updateSelectInput(session, "chr", selected = selected_chrom, choices = visible_chroms)
      shiny::updateSelectizeInput(session, "genes", selected = selected_gene, choices = visible_genes, server = TRUE)
      shiny::updateSelectizeInput(session, "columns", selected = selected_columns, choices = visible_columns)
    }
  )

  # Update list of visible chromosomes when the user changes build
  observeEvent(
    {
      input$build
    },
    {
      shiny::validate(shiny::need(input$password == settings$password, "Password required"))

      visible_chroms <- if (hg_build(input) == "hg19") {
        database$chroms_hg19
      } else {
        database$chroms_hg38
      }

      shiny::updateSelectInput(session, "chr", choices = visible_chroms)

      choices = list()
      choices[[sprintf("Gene (%i)", length(database$genes))]] = "genes"
      choices[[sprintf("Contig (%i)", length(visible_chroms))]] = "contig"

      shiny::updateRadioButtons(session, "select_by", choices = choices, inline = TRUE)
    }
  )

  # This field must NOT be used for anthing but client-side decisions; no data must be
  # sent based on whether or not this is true. It is easily circumvented by the user.
  output$show_ui <- reactive({
    input$password == settings$password
  })
  # Required for the output value to be included despite there being no visible widget
  outputOptions(output, "show_ui", suspendWhenHidden = FALSE)

  # Reset position when the contig is changed or the reset button is pressed
  observeEvent(
    {
      input$chr
      input$btn_reset
    },
    {
      shiny::updateNumericInput(session, "min_pos", value = 1)
      shiny::updateNumericInput(session, "max_pos", value = numeric())
    }
  )

  # Reset filters on request
  observeEvent(
    {
      input$btn_reset
    },
    {
      shiny::updateCheckboxInput(session, "require_pass", value = TRUE)
      shiny::updateSelectInput(session, "consequence", selected = "Any consequence")
      shiny::updateNumericInput(session, "min_maf", value = 0)
      shiny::updateNumericInput(session, "max_maf", value = 1)
      shiny::updateTextAreaInput(session, "query", value = "")
    }
  )

  # Set visible columns on request
  observeEvent(
    {
      input$btn_view_all
    },
    {
      shiny::validate(shiny::need(input$password == settings$password, "Password required"))

      visible_columns <- c("Chr", "Pos", database$columns)
      shiny::updateSelectizeInput(session, "columns", selected = visible_columns, choices = visible_columns)
    }
  )

  # Reset visible columns on request
  observeEvent(
    {
      input$btn_view_std
    },
    {
      shiny::validate(shiny::need(input$password == settings$password, "Password required"))

      visible_columns <- c("Chr", "Pos", database$columns)
      shiny::updateSelectizeInput(session, "columns", selected = settings$columns, choices = visible_columns)
    }
  )

  data <- reactive({
    shiny::validate(shiny::need(input$password == settings$password, "Password required"))

    region <- select_regions(input)
    # Ensure that only valid column names are used
    valid_columns <- c("Chr", "Pos", database$columns)
    visible_columns <- intersect(valid_columns, input$columns)

    if (has_values(region$chr, region$min_pos, visible_columns)) {
      query <- build_query(input, region, params, visible_columns)
      result <- database$query(query$string, params = query$params)

      # Fill in names and setup sorting of consequence terms based on rank
      result <- create_sort_order(result, query$columns, "Func_most_significant")
      result <- create_sort_order(result, query$columns, "Func_least_significant")
      result <- create_sort_order(result, query$columns, "Func_most_significant_canonical")

      result
    }
  })

  observeEvent(
    {
      data
      input$columns
    },
    {
      consequence_columns <- c(
        "Func_most_significant",
        "Func_least_significant",
        "Func_most_significant_canonical"
      )

      # Ensure that only valid column names are used
      valid_columns <- c("Chr", "Pos", database$columns)
      visible_columns <- intersect(valid_columns, input$columns)
      visible_consequences <- intersect(consequence_columns, visible_columns)

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
        {
          results <- data()

          # Remove any columns added for house-keeping
          hidden_columns <- sprintf("%s_order", visible_consequences)
          results <- results[, c(visible_columns, hidden_columns), drop = FALSE]

          results
        },
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
    }
  )

  observeEvent(
    {
      input$password
      input$columns
    },
    {
      shiny::validate(shiny::need(input$password == settings$password, "Password required"))

      names <- c("Chr", "Pos", database$columns_info$Name)
      descs <- c("Contig in the selected build", "Position in the selected build", database$columns_info$Description)

      info <- data.frame(
        Enabled = ifelse(names %in% input$columns, "✓", ""),
        Name = names,
        Description = descs
      )

      output$columns <- DT::renderDataTable(
        info,
        selection = "single",
        options = list(
          scrollX = TRUE,
          scrollY = "80vh",
          paging = FALSE
        )
      )
    }
  )


  output$btn_download <- downloadHandler(
    filename = function() {
      strftime(Sys.time(), "annovep_%Y%m%d_%H%M%S.tsv")
    },
    content = function(file) {
      results <- data()

      # Remove any columns added for house-keeping
      results <- results[, input$columns, drop = FALSE]

      write.table(results, file, quote = FALSE, sep = "\t", row.names = FALSE)
    }
  )

  output$query_errors <- renderUI({
    if (has_values(user_filter$errors)) {
      span(
        HTML("<h5 style='color: #AB0000; text-align: center;'>"),
        as.character(user_filter$errors),
        HTML("</h5>")
      )
    }
  })

  # Constant since this is only for startup errors
  output$database_errors <- reactive({
    !is.null(database$errors)
  })
  # Required for the output value to be included despite there being no visible widget
  outputOptions(output, "database_errors", suspendWhenHidden = FALSE)

  output$db_errors <- renderUI({
    span(
      HTML("<h5 style='color: #AB0000; text-align: center;'>"),
      as.character(database$errors),
      HTML("</h5>")
    )
  })

  observeEvent(
    {
      input$password
      input$table_cell_clicked
    },
    {
        output$json <- renderUI({
        shiny::validate(shiny::need(input$password == settings$password, "Password required"))

        if (length(input$table_cell_clicked) > 0) {
          results <- data()
          row <- results[input$table_cell_clicked$row, , drop=FALSE]

          json <- database$query_vec(
            "SELECT HEX(Data) FROM [JSON] Where [Hg38_chr] = :chr AND [Hg38_pos] = :pos",
            params = list(chr = row$Hg38_chr, pos = row$Hg38_pos)
          )

          json <- paste("[", paste(sprintf('"%s"', json), collapse = ", "), "]")
        } else {
          json <- "null"
        }

        message <- ifelse(
          database$has_json,
          "Click a row to view raw annotation data.",
          "Raw annotation data not available."
        )

        shiny::div(
          shiny::pre(
            shiny::code(class="language-json")
          ),
          # JSON data is stored in compact form and needs to be pretty-printed\n
          shiny::tags$script(HTML(sprintf("
            var input = %s;
            var output = '\"%s\"';

            if (input !== null) {
              var output = [];
              for (var i = 0; i < input.length; ++i) {
                /* Convert hex-encoded data to binary */
                var enc_data = input[i];
                var raw_data = new Uint8Array(enc_data.length / 2);
                for (var j = 0; j < enc_data.length; j += 2) {
                  raw_data[j / 2] = parseInt(enc_data.substr(j, 2), 16);
                }

                /* zlib decompress binary data */
                var data = pako.inflate(raw_data);
                /* Convert binary data to text and parse JSON */
                var text = new TextDecoder().decode(data);
                output.push(JSON.parse(text));
              }

              /* Pretty print JSON */
              output = JSON.stringify(output, null, 2);
            }

            var elem = document.getElementsByClassName(\"language-json\")[0];
            elem.innerHTML = output;

            hljs.highlightAll();
          ", json, message)))
        )
      })
    }
  )

  # Fill in dynamic list of consequence terms
  shiny::updateSelectInput(session, "consequence", choices = c("Any consequence", rev(database$consequences$Name)))

  # Enable automatic reconnections
  session$allowReconnect("force")
}
