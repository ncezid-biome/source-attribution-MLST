#' @title Prediction
#'
#' @param query (character) A CSV file with two rows: a header and values for an MLST profile. The header should only have columns with relevant loci and not even an identifier for the genome.
#' @param model (character) Path to single random forest model RDS file
#' @param ncores (integer, default: 1L) Number of cores to use during random forest fit
#'
#' @return prediction The prediction object from `predict.rfsrc()`
#' @export
#'
#' @examples 
#' See [inst/bootstrapRF.R]
prediction <- function(query, model, ncores = 1L) {

  log_info(paste0("Running with ", ncores, " cores."))
  log_info(paste0("Will read model: ", model))

  orgOpt <- options()
  options(rf.cores = ncores, mc.cores = ncores)
  on.exit(options(orgOpt))

  loci_start_with <- "LMO"

  query <- read.csv(query) %>%
    mutate_all(~ ifelse(is.na(.), 0, .)) %>%
    mutate(across(everything(), ~ as.factor(as.character(.x)))) %>%
    as.data.frame() # rfsrc() doesn't work with a tibble

  m <- readRDS(model)

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

