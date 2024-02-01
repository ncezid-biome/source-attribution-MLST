suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("gt"))
suppressPackageStartupMessages(library("randomForestSRC"))
suppressPackageStartupMessages(library("tidyverse"))

log_threshold(SUCCESS)

test_that("Prediction on LMO0003 with example_query", {
  ncores <- 1

  # rfsrc prediction objects
  predictions <- list()
  for(i in seq(1,3)){
    model <- paste0(rds_dir,"/bs",(i+22),".rds")
    predictions[[i]] <- prediction(model_filename = model, 
                                   query = example_query, ncores = ncores)
    
    
    # Save this file for downstream tests
    out_file <- paste0(rds_dir, "/predictions", (i+22), ".rds")
    #log_success(paste0("Saving file to ", out_file))
    saveRDS(predictions[[i]], file = out_file)
    #log_success("BREAK to quickly dev downstream"); break;
  }

  # Test on the first prediction
  pred <- readRDS(paste0(rds_dir, "/predictions23.rds"))
  my_table <- pred$predicted

  expect_equal(sort(colnames(my_table)), sort(c("dairy", "meat", "vegetable", "fruit", "seafood")), expected.label = "column names on the prediction table")

  # Dairy
  obs <- round(my_table[,"dairy"],3) + 0
  exp <- round(0.1508644,3) + 0
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
  obs <- round(my_table[,"meat"],3) + 0
  exp <- round(0.3109127,3) + 0
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
  obs <- round(my_table[,"vegetable"],3) + 0
  exp <- round(0.2108486,3) + 0
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
  obs <- round(my_table[,"fruit"],3) + 0
  exp <- round(0.1569874,3) + 0
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
  obs <- round(my_table[,"seafood"],3) + 0
  exp <- round(0.1703869,3) + 0
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

