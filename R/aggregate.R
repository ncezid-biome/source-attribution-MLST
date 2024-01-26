#' Aggregate predictions
#' 
#' This function aggregates predictions from multiple
#' calls of `prediction()`.
#'
#' @param predictions (vector of predict.rfsrc() objects) The prediction objects from `prediction()`.
#' 
#' @return something 
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' # Example usage:
#' prediction <- list()
#' predictions[[1]] <- prediction(model = "bs23.rds", ...)
#' predictions[[2]] <- prediction(model = "bs24.rds", ...)
#' result <- aggregate(predictions)
#' }
#' 
#' @importFrom logger log_info
#' @importFrom utils read.csv write.table
#' @importFrom magrittr `%>%`
#' @importFrom randomForestSRC predict.rfsrc

aggregate <- function(predictions) {

  # Get all the confidence scores
  scores <- lapply(predictions, function(pred) pred$predicted)

  # Get averages of different sources
  all_scores <- do.call(rbind, scores)
  
  # Calculate the averages and standard deviations for each category
  category_summary <- apply(all_scores, 2, function(x) c(mean = mean(x), sd = sd(x)))

# TODO return averaged VIMP list
# TODO make the return variable an object holding all this stuff
# TODO document this better

  return(category_summary)
}

