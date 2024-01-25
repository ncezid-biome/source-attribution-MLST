#!/usr/bin/env Rscript

#'bootstrapRF.R
#'
#'Create a random forest model from MLST profiles. The profiles spreadsheet
#'needs to have a column for the genome name, "food", and the MLST profiles on
#'subsequent columns. Run with `--help` for up to date usage.
#'
#'@param input Spreadsheet describing MLST profiles, in csv or csv.gz format.
#'@param output Directory of bootstrap random forest models to output.
#'@param dependent The dependent variable in the spreadsheet. Default:food
#'@param core-loci A comma-separated list of core loci to help remove duplicate 
#'  isolates.
#'  These loci must be present as headers in the spreadsheet from `--input`.
#'@param starts-with The prefix of all independent variables. Default:LMO
#'@param seed Random seed. Default:23
#'@param bootstraps How many random forest bootstrap models to output
#'@param threads How many cores to use. Default:1
#'@param help Show the most up to date usage for this script.
#'
#' @examples
#'
#' # Basic example
#'
#' In this example, we make one random forest model under `results/bs23.rds`,
#' and then we query with the provided example MLST profile
#' `data/example_query.csv`.
#' If more bootstraps were requested than 1, there would be more files under
#' `results/`.
#' Each random forest number is given a bootstrap number in the file name,
#' with the random seed used.
#' Seeds are incremented by 1 if more bootstraps are requested.
#'
#' ```shell
#' Rscript scripts/bootstrapRF.R \
#'   --input data/isolates_original_plus_new_dec_1_2021.csv.gz \
#'   -o results \
#'   --dependent food \
#'   --core-loci data/cgMLST_loci.csv \
#'   --starts-with LMO \
#'   --bootstraps 1 \
#'   --threads 8 \
#'   --seed 23
#' ```
#'
#' # Advanced example
#'
#' You can use `--threads` to use multiple cores in the underlying libraries.
#' However, you can let the operating system parallelize the calls if you
#' have many models to create.
#'
#' Using this exact method, we call 4 random seeds with shuf and parallelize
#' with xargs using 4 processors. The seed is transferred to the program with
#' `--seed`.
#' We request 10 models using `--bootstraps` and therefore we are getting
#' 4 x 10 = 40 models.
#' Formally, it is possible that the seeds are close enough to each other that
#' some output models will override each other.
#' You can try to avoid that by increasing the range in shuf or by creating
#' one-off models to fill in the gaps.
#'
#' ```shell
#' shuf -i 1-9999 -n 4 | \
#'    xargs -n 1 -P 4 bash -c ' \
#'      Rscript scripts/bootstrapRF.R \
#'        --input data/isolates_original_plus_new_dec_1_2021.csv.gz \
#'        -o results/LMO \
#'        --dependent food \
#'        --core-loci data/cgMLST_loci.csv \
#'        --starts-with LMO \
#'        --bootstraps 10 \
#'        --threads 1 \
#'        --seed $0
#'      '
#'``` 
#'
#'@author Lee Katz
#'@keywords MLST random-forest source-attribution food
#'@references <https://doi.org/10.1089/fpd.2023.0046>
#'
  
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logger))

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

sourcerer::bootstrapping(input = opt$input, output = opt$output, core_loci = opt$'core-loci',
              ncores = opt$threads, bootstrap_reps = opt$bootstraps,
              loci_start_with = opt$'starts-with', my_seed = opt$seed )
