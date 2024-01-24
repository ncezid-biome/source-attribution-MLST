#' prediction 
#'
#' @param opt A set of options from optparse
#'
#' @return prediction The prediction object from `predict.rfsrc()`
#' @export
#'
#' @examples 
#' See [inst/bootstrapRF.R]
prediction <- function(opt) {
  ncores <- opt$threads

  model_filename <- opt$model
  log_info(paste0("Running with ",ncores," cores."))
  log_info(paste0("Will read model: ", model_filename))

  orgOpt <- options()
  options(browser = 'firefox') 
  options(rf.cores=ncores,mc.cores=ncores)

  loci_start_with <- "LMO"

  query <- read.csv(opt$query) %>%
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
  pred <- predict.rfsrc(m, newdata = query_filtered_cols)
  return(pred)

  write.table(pred$predicted,
              file = stdout(),
              sep  = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE )

}

