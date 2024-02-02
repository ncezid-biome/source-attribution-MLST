
test_that("Prediction on LMO0003 with example_query", {
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

  # rfsrc prediction objects
  predictions <- list()

  for (i in seq_along(filenames)) {
    model <- filenames[i]
    predictions[[i]] <- prediction(
      model_filename = model,
      query = example_query, ncores = ncores
    )
  }

  # Test on the first prediction
  pred <- predictions[[1]]
  my_table <- pred$predicted

  expect_equal(sort(colnames(my_table)), sort(c("dairy", "meat", "vegetable", "fruit", "seafood")), expected.label = "column names on the prediction table")

  # Dairy
  obs <- round(my_table[, "dairy"], 3) + 0
  exp <- round(0.1508644, 3) + 0
  names(obs) <- "dairy"
  names(exp) <- "dairy"

  expect_equal(
    obs,
    exp,
    info = "Test for dairy percentage",
    label = paste0("Observed dairy percentage (", obs, ")"),
    expected.label = paste0("Expected dairy percentage (", exp, ")")
  )

  # Meat
  obs <- round(my_table[, "meat"], 3) + 0
  exp <- round(0.3109127, 3) + 0
  names(obs) <- "meat"
  names(exp) <- "meat"

  expect_equal(
    obs,
    exp,
    info = "Test for meat percentage",
    label = paste0("Observed meat percentage (", obs, ")"),
    expected.label = paste0("Expected meat percentage (", exp, ")")
  )

  # Vegetable
  obs <- round(my_table[, "vegetable"], 3) + 0
  exp <- round(0.2108486, 3) + 0
  names(obs) <- "vegetable"
  names(exp) <- "vegetable"

  expect_equal(
    obs,
    exp,
    info = "Test for vegetable percentage",
    label = paste0("Observed vegetable percentage (", obs, ")"),
    expected.label = paste0("Expected vegetable percentage (", exp, ")")
  )

  # Fruit
  obs <- round(my_table[, "fruit"], 3) + 0
  exp <- round(0.1569874, 3) + 0
  names(obs) <- "fruit"
  names(exp) <- "fruit"

  expect_equal(
    obs,
    exp,
    info = "Test for fruit percentage",
    label = paste0("Observed fruit percentage (", obs, ")"),
    expected.label = paste0("Expected fruit percentage (", exp, ")")
  )

  # Seafood
  obs <- round(my_table[, "seafood"], 3) + 0
  exp <- round(0.1703869, 3) + 0
  names(obs) <- "seafood"
  names(exp) <- "seafood"

  expect_equal(
    obs,
    exp,
    info = "Test for seafood percentage",
    label = paste0("Observed seafood percentage (", obs, ")"),
    expected.label = paste0("Expected seafood percentage (", exp, ")")
  )
})
