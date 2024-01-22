#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

cli_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", cli_args, value = TRUE)
script_name <- sub("--file=", "", script_arg)

script_name_absolute <- tools::file_path_as_absolute(script_name)
script_dir <- dirname(script_name_absolute)

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

ncores <- opt$threads

model_filename <- opt$model
#print(paste0("Running with ",ncores," cores."))
#print(paste0("Will read model: ", model_filename))

orgOpt <- options()
options(browser = 'firefox') 
options(rf.cores=ncores,mc.cores=ncores)

loci_start_with <- "LMO"

source(paste0(script_dir, "/wgMLST_funs_update.R"))

query <- read.csv(opt$query) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  mutate(across(everything(), ~ as.factor(as.character(.x)))) %>%
  as.data.frame() # rfsrc() doesn't work with a tibble

m <- readRDS(model_filename)

# Identify missing columns
missing_columns <- setdiff(names(query), names(m$xvar))

# Remove missing columns from the query
query_filtered_cols <- query[, setdiff(names(query), missing_columns)]

# Align factor levels
for (col in grep(loci_start_with, names(query_filtered_cols))) {
  query_filtered_cols[[col]] <- factor(
    query_filtered_cols[[col]],
    levels = intersect(levels(query_filtered_cols[[col]]), levels(m$xvar[[col]]))
  )
}

# Make predictions
pred <- predict.rfsrc(m, newdata = query_filtered_cols)

#print(pred$predicted)

#cat(paste(names(pred$predicted), collapse="\t"))
#cat("\n")
#cat(paste(pred$predicted, collapse="\t"))
#cat("\n")

write.table(pred$predicted,
            file = stdout(),
            sep  = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE )

