bootstrapping <- function(opt) {

  bootstrap_folder <- opt$output

  ncores <- opt$threads
  bootstrap_reps <- opt$bootstraps 
  log_info(paste0("Running with ",ncores," cores and ", bootstrap_reps, " bootstraps"))

  orgOpt <- options()
  options(browser = 'firefox') 
  options(rf.cores=ncores,mc.cores=ncores)
  #options(rf.cores=orgOpt$rf.cores,mc.cores=orgOpt$mc.cores)

  loci_start_with <- opt$'starts-with'
  log_info(paste0("Loci start with ",loci_start_with))

  log_info(paste0("Getting the MLST profiles from ",opt$input))
  lm_dat <- read.csv(opt$'input', header = TRUE)
  # This is used in the sel_rep_iso function to select representative isolates
  cgmlst_loci <- read.csv(opt$'core-loci') %>% names

  ### ht defines the threshold of the proportional difference within which isolates were treated
  #### as originated from the same outbreaks or the collection from the same facilities

  ht <- 0.004

  si <- sel_rep_iso(lm_dat, ht, cgmlst_loci) ### select representative isolates


  ###################################################################
  #  importance of genes from random forest model based on all genes
  ###################################################################
  log_info(paste0("Filtering input data to remove nearly duplicate profiles and to just view relevant loci"))
  train.df.all <- lm_dat %>%
    filter(SRR_ID %in% si) %>%
    select("food", starts_with(loci_start_with)) %>%
    mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
    mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
    mutate(food = as.factor(as.character(food))) %>%
    as.data.frame() # rfsrc() doesn't work with a tibble

  # Bootstrapping
  if(!dir.exists(bootstrap_folder)){
    dir.create(bootstrap_folder)
  }
  log_info(paste0("Running bootstraps and saving them to ", bootstrap_folder, "/*.rds"))
  # Set the inital seed for RF models
  my_seed <- opt$seed
  for (i in 1:bootstrap_reps){
    #my_seed <- sample(1:as.integer(.Machine$integer.max))
    log_info(paste0("Modeling rep ", i, " with seed ", my_seed, "..."))
    model <- rfsrc(food ~ ., train.df.all, importance = T, seed = my_seed) #

    # Save intermediate results
    # TODO in the future in might be nice to save the filename with the random seed or a hash
    # so that we can just add more bootstraps if needed
    filename <- paste0(bootstrap_folder,"/bs", my_seed, ".rds")
    log_info(paste0("Saving bootstrap", i, " to ", filename))
    saveRDS(model, file = filename)

    # Free up memory
    #rm(model)
    
    # Next seed: it's not perfect but just increment the seed
    my_seed <- my_seed +1
    set.seed(my_seed)
  }

}

