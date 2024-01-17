# Revised code by Weidong Gu December 8, 2023
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(hierfstat))
suppressPackageStartupMessages(library(poppr))
suppressPackageStartupMessages(library(mmod))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(multiROC))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(msgl))
suppressPackageStartupMessages(library(mlogit))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(StatMatch))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(vcd))
suppressPackageStartupMessages(library(tidyverse))

base_dir <- "/home/gzu2/src/source-attribution-MLST/"
data_folder <- paste0(base_dir,"data")
results_folder <- paste0(base_dir,"results")
plot_folder <- paste0(base_dir,"plots")
script_folder <- paste0(base_dir, "scripts")
#ncores <- ifelse(detectCores() - 1 < 1, 1, detectCores() - 1)
ncores <- 2
bootstrap_reps <- 10
reproducible <- TRUE
print(paste0("Running with ",ncores," cores and ", bootstrap_reps, " replicates. Reproducibility has been set to ", reproducible));

orgOpt <- options()
options(browser = 'firefox') 
options(rf.cores=ncores,mc.cores=ncores)
#options(rf.cores=orgOpt$rf.cores,mc.cores=orgOpt$mc.cores)

loci_start_with <- "LMO"
#print("DEBUG: changed filter from starts_with to starts_with('LMO0030') to reduce it to 10 loci"); loci_start_with <- "LMO0030"
print(paste0("Loci start with ",loci_start_with))

source(paste0(script_folder,"/wgMLST_funs_update.R"))

#lm_dat <- readRDS(paste0(data_folder,"/isolates_original_plus_new_dec_1_2021.rds"))
lm_dat <- read.csv(paste0(data_folder,"/isolates_original_plus_new_dec_1_2021.csv.gz"),header = TRUE)
cgmlst_loci <- read.csv(paste0(data_folder, "/cgMLST_loci.csv")) %>% names

### ht defines the threshold of the proportional difference within which isolates were treated
#### as originated from the same outbreaks or the collection from the same facilities

ht <- 0.004

si <- sel_rep_iso(lm_dat, ht) ### select representative isolates


###################################################################
#  importance of genes from random forest model based on all genes
###################################################################
train.df.all <- lm_dat %>%
  filter(SRR_ID %in% si) %>%
  select("food", starts_with(loci_start_with)) %>%
  mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
  mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
  as.data.frame() # rfsrc() doesn't work with a tibble

# Bootstrapping
model_list <- list()
filenames  <- list()
for (i in 1:bootstrap_reps){
  model <- rfsrc(food ~ ., train.df.all, importance = T ) #

  # Save intermediate results
  filename <- paste0(results_folder,"/model", i, ".rds")
  print(paste0("Saving bootstrap", i, " to ", filename))
  saveRDS(model, file = filename)

  # Remember the intermediate files
  model_list[[i]] <- readRDS(filename)
  filenames[[i]] <- filename
}

# If we made it this far, then combine all the models into an array of models
saveRDS(model_list, paste0(results_folder, "/models.rds"))
# Delete intermediate files
for (f in filenames) {
  print(paste0("Removing intermediate file ", f))
  file.remove(f)
}

#rf.all <- rfsrc(food ~ ., train.df.all, importance = T ) #

#rf <- rf.all
#c.imp <- rf$importance[, 1]
#rk.genes <- names(sort(c.imp, decreasing = T)) # [ZC: top ranked genes]
#c.imp <- as.numeric(c.imp)

#for (i in seq_along(rk.genes)) {
#    cat(rk.genes[i], "\t", c.imp[i], "\n")
#}



