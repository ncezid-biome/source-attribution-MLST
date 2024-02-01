#' @title Listeria isolates
#' 
#' @description This dataset provides an example spreadsheet for creating 
#' a random forest model.
#' 
#' @name Listeria_isolates
#' 
#' @export Listeria_isolates
#' 
#' @format A file path to Listeria_isolates.csv.gz
Listeria_isolates <- system.file("data", "Listeria_isolates.csv.gz", package = "sourcerer")

#' @title Example query
#' 
#' @description This is an example query for use with the prediction step.
#' 
#' @name example_query
#' 
#' @export example_query
#' 
#' @format A file path to example_query.csv
example_query <- system.file("data", "example_query.csv", package = "sourcerer")

#' @title rds folder
#' 
#' @description This is where RDS files are stored for the unit tests
#' 
#' @name rds_dir
#' 
#' @export rds_dir
#' 
#' @format A file path to the rds_dir
rds_dir <- system.file("tests", "test-results", package = "sourcerer")
