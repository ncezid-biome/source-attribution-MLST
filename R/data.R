#' @title Listeria isolates
#' @description This dataset provides an example spreadsheet for creating
#' a random forest model.
#' @format A data frame with 719 rows and 4807 variables:
#' \describe{
#'   \item{\code{SRR_ID}}{character Sample ID}
#'   \item{\code{SourceSite}}{character Sample Source}
#'   \item{\code{food}}{character Food Type}
#'   \item{\code{LMO0XXX}{double Expression at site XXX}}
#' }
"Listeria_isolates"


#' @title Example query
#' @description This is an example query for use with the prediction step.
#' @format A data frame with 1 rows and 4804 variables:
#' \describe{
#'   \item{\code{LMO0XXX}}{double Expression at site XXX}
#' }
"example_query"
