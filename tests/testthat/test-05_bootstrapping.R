suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("gt"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("randomForestSRC"))
suppressPackageStartupMessages(library("tidyverse"))



test_that("Bootstrapping on LMO0003", {

  # Actually run the options but accept some slightly modified defaults.
  # Therefore, all options have been set without user input.
  option_list <- list(
      optparse::make_option(c("-i", "--input"), type = "character", help = "Spreadsheet describing MLST profiles, in csv or csv.gz format.", default = "Listeria_isolates.csv.gz"),
      optparse::make_option(c("-o", "--output"), type = "character", help = "Directory of bootstrap random forest models to output", default = "test-results"),
      optparse::make_option(c("-d", "--dependent"), type = "character", help = "The dependent variable in the spreadsheet. Default:food", default = "food"),
#      optparse::make_option(c("-c", "--core-loci"), type = "character", help = "A comma-separated list of core loci to help remove duplicate isolates. These loci must be present as headers in the spreadsheet from --input.", default = "cgMLST_loci.csv"),
      optparse::make_option(c("", "--starts-with"), type = "character", help = "The prefix of all independent variables. Default:LMO", default = "LMO0003"),
      optparse::make_option(c("", "--seed"), type = "integer", help = "Random seed. Default:23", default = 23),
      optparse::make_option(c("-b", "--bootstraps"), type = "integer", help = "How many random forest bootstraps to output", default = 3),
      optparse::make_option(c("-t", "--threads"), type = "integer", help = "How many cores to use. Default:1", default = 1)
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)

  filenames <- bootstrapping(input = opt$input, output = opt$output, 
                  ncores = opt$threads, bootstrap_reps = opt$bootstraps,
                  loci_start_with = opt$'starts-with', my_seed = opt$seed )

  #filenames <- bootstrapping(opt)

  expect_equal(filenames[[1]], "test-results/bs23.rds", info = "Filename 1")
  expect_equal(filenames[[2]], "test-results/bs24.rds", info = "Filename 2")
  expect_equal(filenames[[3]], "test-results/bs25.rds", info = "Filename 3")


  rf      <- list()
  rf[[1]] <- readRDS(filenames[[1]])
  expect_equal(rf[[1]]$forest$n, 719, info = "Number of alleles in the first RF model")
  expect_equal(rf[[1]]$forest$yvar.names, "food", info = "yvar in the first RF model")
  expect_equal(rf[[1]]$importance["LMO00030", "all"], 0.07586199, tolerance = 0.01, info = "importance of LMO00030 in the all category")

})

