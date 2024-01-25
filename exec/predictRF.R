#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logger))

# Command line argument parsing
option_list <- list(
    make_option(c("-m", "--model"), type = "character", help = "A single random forest model RDS file"),
    make_option(c("-q", "--query"), type = "character", help = "A CSV file with two rows: a header and values for an MLST profile. The header should only have columns with relevant loci and not even an identifier for the genome."),
    make_option(c("-t", "--threads"), type = "integer", help = "How many cores to use. Default: 1", default = 1)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# required options
required_options <- c("model", "query")
for (o in required_options) {
  if (!(o %in% names(opt))) {
    cat("ERROR: Required option", o, "is missing.\n")
    print_help(opt_parser)
    q(status = 1)
  }
}

sourcerer::prediction(model_filename = opt$model, 
                  query = opt$query, ncores = opt$threads)
