# Revised code by Weidong Gu December 8, 2023
suppressPackageStartupMessages(library(cluster))
#suppressPackageStartupMessages(library(hierfstat))
#suppressPackageStartupMessages(library(poppr))
#suppressPackageStartupMessages(library(mmod))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(doParallel))
#suppressPackageStartupMessages(library(multiROC))
suppressPackageStartupMessages(library(gt))
#suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(geiger))
#suppressPackageStartupMessages(library(dendextend))
#suppressPackageStartupMessages(library(DECIPHER))
#suppressPackageStartupMessages(library(msgl))
#suppressPackageStartupMessages(library(mlogit))
#suppressPackageStartupMessages(library(caret))
#suppressPackageStartupMessages(library(xgboost))
#suppressPackageStartupMessages(library(Matrix))
#suppressPackageStartupMessages(library(StatMatch))
#suppressPackageStartupMessages(library(grid))
#suppressPackageStartupMessages(library(gtable))
#suppressPackageStartupMessages(library(vcd))
suppressPackageStartupMessages(library(tidyverse))

base_dir <- "/home/gzu2/src/source-attribution-MLST/"
data_folder <- paste0(base_dir,"data")
results_folder <- paste0(base_dir,"results")
plot_folder <- paste0(base_dir,"plots")
script_folder <- paste0(base_dir, "scripts")
#ncores <- ifelse(detectCores() - 1 < 1, 1, detectCores() - 1)
ncores <- 1
reproducible <- TRUE
model_filename <- paste0(results_folder,"/models.rds")
print(paste0("Running with ",ncores," cores. Reproducibility has been set to ", reproducible));
print(paste0("Will read model with bootstraps: ", model_filename))

orgOpt <- options()
options(browser = 'firefox') 
options(rf.cores=ncores,mc.cores=ncores)

loci_start_with <- "LMO"
#print("DEBUG: changed filter from starts_with to starts_with('LMO0030') to reduce it to 10 loci"); loci_start_with <- "LMO0030"

source(paste0(script_folder,"/wgMLST_funs_update.R"))

lm_dat <- read.csv(paste0(data_folder,"/isolates_original_plus_new_dec_1_2021.csv.gz"),header = TRUE)
cgmlst_loci <- read.csv(paste0(data_folder, "/cgMLST_loci.csv")) %>% names

### ht defines the threshold of the proportional difference within which isolates were treated
#### as originated from the same outbreaks or the collection from the same facilities

ht <- 0.004

si <- sel_rep_iso(lm_dat, ht) ### select representative isolates

query <- lm_dat %>%
  select(starts_with(loci_start_with)) %>%
  mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
  mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
  as.data.frame() # rfsrc() doesn't work with a tibble
query <- query[1,]; print("DEBUG: using first row of spreadsheet as query"); 

model <- readRDS(model_filename)
print(paste0("Read models file ", model_filename, " which has ", length(model), " replicates"))

predictions <- list()
for (i in 1:length(model)) {

  # Save the model
  m <- model[[i]]

	missing_columns <- setdiff(names(query), names(m$xvar))
	print(missing_columns)

	query_filtered_cols <- query[, setdiff(names(query), missing_columns)]

	# Sort factor levels in the training data
	m$xvar[, grep(loci_start_with, names(m$xvar))] <- 
		lapply(m$xvar[, grep(loci_start_with, names(m$xvar))], function(x) factor(x, levels = sort(levels(x))))
	# Sort factor levels in the new data
	query_filtered_cols[, grep(loci_start_with, names(query_filtered_cols))] <- 
		lapply(query_filtered_cols[, grep(loci_start_with, names(query_filtered_cols))], function(x) factor(x, levels = sort(levels(x))))

	# Check factor levels in the training data
	levels_in_training <- lapply(m$xvar[, grep(loci_start_with, names(m$xvar))], levels)
	# Check factor levels in the new data
	levels_in_new_data <- lapply(query_filtered_cols[, grep(loci_start_with, names(query_filtered_cols))], levels)
	# Compare factor levels
	print(levels_in_training)
	print(levels_in_new_data)

	# Identify extra levels in the test data
	extra_levels_test <- lapply(
		query_filtered_cols[, grep(loci_start_with, names(query_filtered_cols))],
		function(x) setdiff(levels(x), unlist(levels_in_training))
	)

	# Remove extra levels in the test data
	for (col in grep(loci_start_with, names(query_filtered_cols))) {
		query_filtered_cols[[col]] <- factor(
			query_filtered_cols[[col]],
			levels = intersect(levels(query_filtered_cols[[col]]), levels_in_training[[col]])
		)
	}

	pred <- predict.rfsrc(m, newdata = query_filtered_cols)
  predictions[[i]] <- pred$predicted
	#print(pred)
  #
	#print("SUMMARY OF PREDICTIONS")
	#print(summary(pred))
	#print("STR OF PREDICTIONS")
	#print(str(pred))
	#print("PREDICTIONS$predicted")
	#print(pred$predicted)
	#print(summary(pred$predicted))
	#print(str(pred$predicted))

	#all_pred <- lapply(model, function(m) predict(m, newdata = query))
	#print(all_pred)
}

print("PREDICTIONS ARRAY")
print(predictions)

