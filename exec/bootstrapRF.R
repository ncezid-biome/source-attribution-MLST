#!/usr/bin/env Rscript
  
#suppressPackageStartupMessages(library(cluster))
#suppressPackageStartupMessages(library(randomForestSRC))
#suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(dendextend))
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(optparse))
#suppressPackageStartupMessages(library(logger))

devtools::load_all()

# Command line argument parsing
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Spreadsheet describing MLST profiles, in csv or csv.gz format."),
    make_option(c("-o", "--output"), type = "character", help = "Directory of bootstrap random forest models to output", default = "results"),
    make_option(c("-d", "--dependent"), type = "character", help = "The dependent variable in the spreadsheet. Default:food", default = "food"),
    make_option(c("-c", "--core-loci"), type = "character", help = "A comma-separated list of core loci to help remove duplicate isolates. These loci must be present as headers in the spreadsheet from --input."),
    make_option(c("", "--starts-with"), type = "character", help = "The prefix of all independent variables. Default:LMO", default = "LMO"),
    make_option(c("", "--seed"), type = "integer", help = "Random seed. Default:23", default = 23),
    make_option(c("-b", "--bootstraps"), type = "integer", help = "How many random forest bootstraps to output", default = 1),
    make_option(c("-t", "--threads"), type = "integer", help = "How many cores to use. Default:1", default = 1)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# required options
required_options <- c("input", "output", "core-loci")
for (o in required_options) {
  if (!(o %in% names(opt))) {
    cat("ERROR: Required option", o, "is missing.\n")
    print_help(opt_parser)
    q(status = 1)
  }
}

bootstrapping(input = opt$input, output = opt$output, core_loci = opt$'core-loci',
              ncores = opt$threads, bootstrap_reps = opt$bootstraps,
              loci_start_with = opt$'starts-with', my_seed = opt$seed )

