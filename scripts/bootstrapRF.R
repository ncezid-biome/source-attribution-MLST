#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(randomForestSRC))
#suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# Get the script name from command-line arguments
cli_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", cli_args, value = TRUE)
script_name <- sub("--file=", "", script_arg)

script_name_absolute <- tools::file_path_as_absolute(script_name)
script_dir <- dirname(script_name_absolute)

# Command line argument parsing
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Spreadsheet describing MLST profiles, in csv or csv.gz format."),
    make_option(c("-o", "--output"), type = "character", help = "Directory of bootstrap random forest models to output"),
    make_option(c("-d", "--dependent"), type = "character", help = "The dependent variable in the spreadsheet", default = "food"),
    make_option(c("-c", "--core-loci"), type = "character", help = "A comma-separated list of core loci to help remove duplicate isolates. These loci must be present as headers in the spreadsheet from --input."),
    make_option(c("-s", "--starts-with"), type = "character", help = "The prefix of all independent variables.", default = "LMO"),
    make_option(c("-b", "--bootstraps"), type = "integer", help = "How many random forest bootstraps to output", default = 1),
    make_option(c("-t", "--threads"), type = "integer", help = "How many cores to use", default = 1)
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

bootstrap_folder <- opt$output

ncores <- opt$threads
bootstrap_reps <- opt$bootstraps 
print(paste0("Running with ",ncores," cores and ", bootstrap_reps, " bootstraps"))

orgOpt <- options()
options(browser = 'firefox') 
options(rf.cores=ncores,mc.cores=ncores)
#options(rf.cores=orgOpt$rf.cores,mc.cores=orgOpt$mc.cores)

loci_start_with <- opt$'starts-with'
print(paste0("Loci start with ",loci_start_with))

functions_script <- file.path(script_dir, "wgMLST_funs_update.R")
source(functions_script)

print(paste0("Getting the MLST profiles from ",opt$input))
lm_dat <- read.csv(opt$'input', header = TRUE)
# This is used in the sel_rep_iso function to select representative isolates
cgmlst_loci <- read.csv(opt$'core-loci') %>% names

### ht defines the threshold of the proportional difference within which isolates were treated
#### as originated from the same outbreaks or the collection from the same facilities

ht <- 0.004

si <- sel_rep_iso(lm_dat, ht) ### select representative isolates


###################################################################
#  importance of genes from random forest model based on all genes
###################################################################
print(paste0("Filtering input data to remove nearly duplicate profiles and to just view relevant loci"))
train.df.all <- lm_dat %>%
  filter(SRR_ID %in% si) %>%
  select("food", starts_with(loci_start_with)) %>%
  mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
  mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
  as.data.frame() # rfsrc() doesn't work with a tibble

# Bootstrapping
if(!dir.exists(bootstrap_folder)){
  dir.create(bootstrap_folder)
}
print(paste0("Running bootstraps and saving them to ", bootstrap_folder, "/*.rds"))
for (i in 1:bootstrap_reps){
  print(paste0("Modeling rep ", i, "..."))
  model <- rfsrc(food ~ ., train.df.all, importance = T ) #

  # Save intermediate results
  # TODO in the future in might be nice to save the filename with the random seed or a hash
  # so that we can just add more bootstraps if needed
  filename <- paste0(bootstrap_folder,"/bs", i, ".rds")
  print(paste0("Saving bootstrap", i, " to ", filename))
  saveRDS(model, file = filename)

  # Free up memory
  rm(model)
}

