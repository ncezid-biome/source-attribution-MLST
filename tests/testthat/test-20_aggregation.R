suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("gt"))
suppressPackageStartupMessages(library("randomForestSRC"))
suppressPackageStartupMessages(library("tidyverse"))


test_that("Aggregating LMO0003 with example_query", {

  ncores <- 1

  # rfsrc prediction objects
  predictions <- list()
  for(i in seq(1,3)){
    model <- paste0("test-results/bs",(i+22),".rds")
    predictions[[i]] <- prediction(model_filename = model, 
                                   query = "example_query.csv", ncores = ncores)

  }

  bootstrapped_prediction <- aggregate(predictions = predictions)
 
  expected <- matrix(
    c(
      0.1501267932, 0.1578980673, 0.3111404334, 0.1697799469, 0.211054759,
      0.0006860643, 0.0008925544, 0.0002023788, 0.0009216065, 0.000183757
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("mean", "sd"), c("dairy", "fruit", "meat", "seafood", "vegetable"))
  )

  expect_equal(bootstrapped_prediction["mean","dairy"], 
               expected["mean","dairy"], 
               tolerance = 0.1)
  expect_equal(bootstrapped_prediction, expected, tolerance = 0.1)

})
