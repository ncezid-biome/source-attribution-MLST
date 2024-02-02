
test_that("Aggregating model LMO0003", {
  input <- system.file("extdata/Listeria_isolates.csv.gz", package = "sourcerer")
  output <- tempdir()
  ncores <- 1
  bootstrap_reps <- 3
  loci_start_with <- "LMO0003"
  my_seed <- 23

  filenames <- bootstrapping(
    input = input, output = output,
    ncores = ncores, bootstrap_reps = bootstrap_reps,
    loci_start_with = loci_start_with, my_seed = my_seed
  )

  models <- lapply(filenames, readRDS)

  composite_model <- aggregate_model(models)

  expect_s3_class(composite_model, "rfsrc")
})

test_that("Aggregating LMO0003 with example_query", {
  input <- system.file("extdata/Listeria_isolates.csv.gz", pacakge = "sourcerer")
  output <- tempdir()
  ncores <- 1
  bootstrap_reps <- 3
  loci_start_with <- "LMO0003"
  my_seed <- 23

  filenames <- bootstrapping(
    input = input, output = output,
    ncores = ncores, bootstrap_reps = bootstrap_reps,
    loci_start_with = loci_start_with, my_seed = my_seed
  )

  # rfsrc prediction objects
  predictions <- list()

  for (i in seq_along(filenames)) {
    model <- filenames[i]
    predictions[[i]] <- prediction(
      model_filename = model,
      query = example_query, ncores = ncores
    )
  }

  bootstrapped_prediction <- aggregate_predictions(predictions = predictions)

  expected <- matrix(
    c(
      0.1501267932, 0.1578980673, 0.3111404334, 0.1697799469, 0.211054759,
      0.0006860643, 0.0008925544, 0.0002023788, 0.0009216065, 0.000183757
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("mean", "sd"), c("dairy", "fruit", "meat", "seafood", "vegetable"))
  )

  expect_equal(bootstrapped_prediction["mean", "dairy"],
    expected["mean", "dairy"],
    tolerance = 0.1
  )
  expect_equal(bootstrapped_prediction, expected, tolerance = 0.1)
})
