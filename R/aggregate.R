#' Aggregate predictions
#' 
#' This function aggregates predictions from multiple
#' calls of `prediction()`.
#'
#' @param predictions (vector of predict.rfsrc() objects) The prediction objects from `prediction()`.
#' 
#' @return category_summary (data frame) A table of means and standard deviations of confidences for each category
#' 
#' @export
#'
#' @examples 
#' # Example usage:
#' prediction <- list()
#' predictions[[1]] <- prediction(model = "bs23.rds", ...)
#' predictions[[2]] <- prediction(model = "bs24.rds", ...)
#' result <- aggregate_predictions(predictions)
#' print(result)
#' 
#' @importFrom logger log_info
#' @importFrom utils read.csv write.table
#' @importFrom magrittr `%>%`
#' @importFrom randomForestSRC predict.rfsrc

aggregate_predictions <- function(predictions) {

  # Get all the confidence scores
  scores <- lapply(predictions, function(pred) pred$predicted)

  # Get averages of different sources
  all_scores <- do.call(rbind, scores)
  
  # Calculate the averages and standard deviations for each category
  category_summary <- apply(all_scores, 2, function(x) c(mean = mean(x), sd = sd(x)))

# TODO return averaged VIMP list - oops this should be in the model and not prediction

  return(category_summary)
}

#' Aggregate models
#' 
#' This function makes a composite model from many models
#' coming from `bootstrapping()`
#' 
#' @param models (vector of `rfsrc`` objects) 
#' 
#' @return composite (`rfsrc` object) The first rfsrc object is returned
#' but modified with an attribute `$aggregate_rank`:
#' each gene is given the median rank in order of importance
#' according to the VIMP scores from each bootstrap.
#' 
#' @export 
#' 
#' @examples 
#'   models = c(
#'     paste0(rds_dir,"/bs23.rds"),
#'     paste0(rds_dir,"/bs24.rds"),
#'     paste0(rds_dir,"/bs25.rds")
#'   )
#'   ref_model <- aggregate_model(models)
#'   print(ref_model$aggregate_rank)
#'  
#' @importFrom logger log_info
#' @importFrom randomForestSRC rfsrc
#' @importFrom rlang duplicate
#' @importFrom stats setNames
aggregate_model <- function(models) {
  # Deep clone the first model so that we don't mess with it
  # by reference accidentally.
  ref_model <- duplicate(models[[1]], shallow = FALSE)
  
  rank_vec <- vector("list", length = length(models))

  for (i in seq_along(models)) {
    rank_vec[[i]] <- setNames(
      order(models[[i]][["importance"]][,"all"]),
      names(models[[i]][["importance"]][,"all"])
    )
  }
  
  # Remove the first row since we seeded it with the 
  # reference.
  vimp_rank_df <- bind_rows(rank_vec)

  ref_model$aggregate_rank <- vimp_rank_df 

  return(ref_model)
}
