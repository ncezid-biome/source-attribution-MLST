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

test_that("Bootstrapping on LMO0003", {

  input  <- Listeria_isolates 
  output <- rds_dir
  ncores <- 1
  bootstrap_reps <- 3
  loci_start_with <- "LMO0003"
  my_seed <- 23
  filenames <- bootstrapping(input = input, output = output, 
                  ncores = ncores, bootstrap_reps = bootstrap_reps,
                  loci_start_with = loci_start_with, my_seed = my_seed )

  #filenames <- bootstrapping(opt)

  expect_equal(basename(filenames[[1]]), 
      basename("test-results/bs23.rds"), info = "Filename 1")
  expect_equal(basename(filenames[[2]]), 
      basename("test-results/bs24.rds"), info = "Filename 2")
  expect_equal(basename(filenames[[3]]), 
      basename("test-results/bs25.rds"), info = "Filename 3")


  rf      <- list()
  rf[[1]] <- readRDS(filenames[[1]])
  expect_equal(rf[[1]]$forest$n, 719, info = "Number of alleles in the first RF model")
  expect_equal(rf[[1]]$forest$yvar.names, "food", info = "yvar in the first RF model")
  expect_equal(rf[[1]]$importance["LMO00030", "all"], 0.07586199, tolerance = 0.01, info = "importance of LMO00030 in the all category")

})

