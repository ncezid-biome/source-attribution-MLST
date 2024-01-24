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
      optparse::make_option(c("-i", "--input"), type = "character", help = "Spreadsheet describing MLST profiles, in csv or csv.gz format.", default = "../../data/isolates_original_plus_new_dec_1_2021.csv.gz"),
      optparse::make_option(c("-o", "--output"), type = "character", help = "Directory of bootstrap random forest models to output", default = "test-results"),
      optparse::make_option(c("-d", "--dependent"), type = "character", help = "The dependent variable in the spreadsheet. Default:food", default = "food"),
      optparse::make_option(c("-c", "--core-loci"), type = "character", help = "A comma-separated list of core loci to help remove duplicate isolates. These loci must be present as headers in the spreadsheet from --input.", default = "../../data/cgMLST_loci.csv"),
      optparse::make_option(c("", "--starts-with"), type = "character", help = "The prefix of all independent variables. Default:LMO", default = "LMO0003"),
      optparse::make_option(c("", "--seed"), type = "integer", help = "Random seed. Default:23", default = 23),
      optparse::make_option(c("-b", "--bootstraps"), type = "integer", help = "How many random forest bootstraps to output", default = 3),
      optparse::make_option(c("-t", "--threads"), type = "integer", help = "How many cores to use. Default:1", default = 1)
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)


  file_stat <- bootstrapping(opt)
  file_stat <- list()
  for(filename in c("test-results/bs23.rds", "test-results/bs24.rds", "test-results/bs25.rds")){
    file_stat <- c(file_stat, list(file.info(filename)))
  }

  # Assuming file_stat is the list obtained from your previous code
  file_sizes <- sapply(file_stat, function(info) info$size)

  # Expected file sizes
  expected_sizes <- c(5052380, 5055516, 5047648)

  # Test file sizes using expect_equal
  for (i in seq_along(file_sizes)) {
    obs <- as.numeric(file_sizes[i])
    exp <- as.numeric(expected_sizes[i])

    expect_equal(obs, exp,
                 tolerance = 0.001,
                 #info = paste("File:", names(file_sizes)[i]),
                 label = "File Size")
  }

})

