#!/usr/bin/env Rscript
requireNamespace("shiny")

main <- function(args) {
  if (length(args) != 1) {
    cat("Usage:\n")
    cat("Rscript <sqlite3.db>\n")
    return(1)
  }

  Sys.setenv(ANNOVEP_DATABASE = args[1])

  cat("starting server ...\n")

  return(shiny::runApp(shiny::shinyAppDir(".", options = list(port = 56789))))
}


quit(status = main(commandArgs(TRUE)))
