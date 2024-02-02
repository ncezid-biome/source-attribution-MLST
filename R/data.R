#' @title Listeria isolates
#' 
#' @description This dataset provides an example spreadsheet for creating 
#' a random forest model.
#' 
#' @name Listeria_isolates
#' 
#' @export Listeria_isolates
#' 
#' @format A file path to Listeria_isolates.csv.gz.
#' This file contains MLST profiles for many Listeria genomes
#' and their food vehicle.
#' It is used for training the random forest model.
Listeria_isolates <- system.file("extdata", "Listeria_isolates.csv.gz", package = "sourcerer")

#' @title Example query
#' 
#' @description This is an example query for use with the prediction step.
#' It contains one genome and its MLST profile.
#' 
#' @name example_query
#' 
#' @export example_query
#' 
#' @format A file path to example_query.csv
example_query <- system.file("extdata", "example_query.csv", package = "sourcerer")

#' @title rds folder
#' 
#' @description This is where RDS files are stored for the unit tests.
#' This includes but might not be limited to:
#' * random forest bootstrap models
#' * predictions files from `example_query`
#' 
#' @name rds_dir
#' 
#' @export rds_dir
#' 
#' @format A file path to the rds_dir
rds_dir <- dirname(Listeria_isolates)
