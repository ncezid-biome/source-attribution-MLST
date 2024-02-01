#' @title Listeria isolates
#' 
#' @description This dataset provides an example spreadsheet for creating 
#' a random forest model.
#' 
#' @name Listeria_isolates
#' 
#' 
#' @format A file path to Listeria_isolates.csv.gz
#' @importFrom usethis use_data
Listeria_isolates <- system.file("data/Listeria_isolates.csv.gz", package = "sourcerer")
use_data(Listeria_isolates, internal = FALSE, overwrite = TRUE, 
          compress = "gzip", version = 3)

#' @title Example query
#' 
#' @description This is an example query for use with the prediction step.
#' 
#' @name example_query
#' 
#' @format A file path to example_query.csv
#' @importFrom usethis use_data
example_query <- system.file("data/example_query.csv", package = "sourcerer")
use_data(example_query, internal = FALSE, overwrite = TRUE, 
          compress = "gzip", version = 3)
