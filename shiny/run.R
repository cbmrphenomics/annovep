#!/usr/bin/env Rscript
requireNamespace("shiny")

main <- function(args) {
  if (length(args) < 1 || length(args) > 3) {
    cat("Usage:\n")
    cat("Rscript run.R <appDir> [<host>] [<port>]\n")
    return(1)
  }

  app_dir <- args[1]
  host <- ifelse(length(args) > 1, args[2], "0.0.0.0")
  port <- ifelse(length(args) > 2, as.numeric(args[3]), 3838)

  cat("starting server in '", app_dir, "', binding to ", host, ":", port, "\n", sep = "")

  return(shiny::runApp(app_dir, port = port, host = host))
}


quit(status = main(commandArgs(TRUE)))
