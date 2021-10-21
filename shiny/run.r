#!/usr/bin/env Rscript

main <- function(args) {
  if (length(args) != 1) {
    cat("Usage:\n")
    cat("Rscript <sqlite3.db>\n")
    return(1)
  }

  library(shiny)

  Sys.setenv(ANNOVEP_DATABASE = args[1])

  cat("starting server ...\n")

  return(runApp(shinyAppDir(".", options = list(port = 56789))))
}


quit(status = main(commandArgs(TRUE)))
