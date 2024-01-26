#' Predict using a random forest model
#' 
#' This function makes predictions using a random forest model.
#' In many situations, you would want to make predictions from
#' many random forest models and aggregate them downstream. 
#'
#' @param model_filename (character) The filename of the Random Forest model.
#' @param query (character) The filename of the query data in CSV format.
#'   This is an MLST profiles spreadsheet.
#' @param ncores (integer, default: 1L) The number of cores to use for parallel processing.
#'
#' @return prediction The prediction object from `predict.rfsrc()`
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' # Example usage:
#' result <- prediction(model_filename = "results/bs23.rds",
#'                      query = "tests/testthat/example_query.csv", ncores = 4)
#' }
#' 
#' @importFrom utils read.csv write.table
#' @importFrom magrittr `%>%`
#' @importFrom randomForestSRC predict.rfsrc
prediction <- function(model_filename, query, ncores = 1L) {
  log_info(paste0("Running with ",ncores," cores."))
  log_info(paste0("Will read model: ", model_filename))

  orgOpt <- options()
  options(rf.cores = ncores, mc.cores = ncores)
  on.exit(options(orgOpt))

  loci_start_with <- "LMO"

  query <- read.csv(query) %>%
    mutate_all(~ ifelse(is.na(.), 0, .)) %>%
    mutate(across(everything(), ~ as.factor(as.character(.x)))) %>%
    as.data.frame() # rfsrc() doesn't work with a tibble

  m <- readRDS(model_filename)

  # Identify missing columns
  missing_columns <- setdiff(names(query), names(m$xvar))

  # Remove missing columns from the query
  query_filtered_cols <- query[, setdiff(names(query), missing_columns)]

  # Align factor levels
  for (col in grep(loci_start_with, names(query_filtered_cols))) {
    query_filtered_cols[[col]] <- factor(
      query_filtered_cols[[col]],
      levels = intersect(levels(query_filtered_cols[[col]]), levels(m$xvar[[col]]))
    )
  }

  # Make predictions
  pred <- randomForestSRC::predict.rfsrc(m, newdata = query_filtered_cols)

  write.table(pred$predicted,
              file = stdout(),
              sep  = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)

  return(pred)
}

