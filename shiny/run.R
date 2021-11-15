#!/usr/bin/env Rscript
requireNamespace("shiny")

main <- function(args) {
  if (length(args) < 1 || length(args) > 3) {
    cat("Usage:\n")
    cat("Rscript run.R [<appDir>] [<host>] [<port>]\n")
    return(1)
  }

  app_dir <- ifelse(length(args) > 0, args[1], ".")
  host <- ifelse(length(args) > 1, args[2], "127.0.0.1")
  port <- ifelse(length(args) > 2, as.numeric(args[3]), 8000)

  cat("starting server in '", app_dir, "', binding to ", host, ":", port, "\n", sep = "")

  return(shiny::runApp(shiny::shinyAppDir(app_dir, options = list(port = port, host = host))))
}


quit(status = main(commandArgs(TRUE)))
