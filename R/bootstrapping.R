#' Bootstrap random forest analysis
#' 
#' This function creates one or more random forest models from a spreadsheet of MLST profiles.
#'
#' @param input A character string specifying the input file path for MLST profile data in csv or csv.gz format.
#' @param output A character string specifying the output directory for random forest models.
#' @param core_loci A csv or csv.gz file with a comma-separated list of core loci to help remove duplicate isolates.
#'  These loci must be present as headers in the spreadsheet from `--input`.
#' @param ncores An integer specifying the number of cores to use for parallel processing. Default is 1.
#' @param bootstrap_reps An integer specifying the number of bootstrap replicates. Default is 1.
#' @param loci_start_with A character string specifying the prefix for loci names. Default is "LMO".
#' @param my_seed An integer specifying the seed for reproducibility. Default is 23.
#' @param ... Additional arguments for future expansion. Currently not used.
#'
#' @return my_filenames A list of random forest filenames
#' 
#' @export
#'
#' @examples TODO
#' \dontrun{
#' # Example usage:
#' model_filenames <- bootstrapping(input = "tests/testthat/isolates_original_plus_new_dec_1_2021.csv.gz", output = "results/",
#'                          core_loci = "tests/testthat/cgMLST_loci.csv", ncores = 4, bootstrap_reps = 100)
#' }
#'
#' @import cluster
#' @import dplyr
#' @importFrom stats cutree hclust quantile rmultinom sd xtabs
#' @importFrom utils read.csv write.table
#' @importFrom magrittr `%>%`
#' @import ggplot2
bootstrapping <- function(input, output, core_loci, ncores = 1,
                          bootstrap_reps = 1, loci_start_with = "LMO",
                          my_seed = 23) {


  if (missing(input)) {
    stop("input is a mandatory argument")
  }
  if (missing(output)){
    stop("output is a mandatory argument")
  }
  if (missing(core_loci)){
    stop("core_loci is a mandatory argument")
  }

  #ncores <- opt$threads
  #bootstrap_reps <- opt$bootstraps 
  log_info(paste0("Running with ",ncores," cores and ", bootstrap_reps, " bootstraps"))

  orgOpt <- options()
  options(browser = 'firefox') 
  options(rf.cores=ncores,mc.cores=ncores)
  #options(rf.cores=orgOpt$rf.cores,mc.cores=orgOpt$mc.cores)

  #loci_start_with <- opt$'starts-with'
  log_info(paste0("Loci start with ", loci_start_with))

  log_info(paste0("Getting the MLST profiles from ", input))
  lm_dat <- read.csv(input, header = TRUE)
  # This is used in the sel_rep_iso function to select representative isolates
  cgmlst_loci <- read.csv(core_loci) %>% names

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
  if(!dir.exists(output)){
    dir.create(output)
  }
  log_info(paste0("Running bootstraps and saving them to ", output, "/*.rds"))
  # Set the inital seed for RF models
  #my_seed <- opt$seed
  my_filenames <- list()
  for (i in 1:bootstrap_reps){
    #my_seed <- sample(1:as.integer(.Machine$integer.max))
    log_info(paste0("Modeling rep ", i, " with seed ", my_seed, "..."))
    model <- rfsrc(food ~ ., train.df.all, importance = T, seed = my_seed) #

    # Save intermediate results
    # TODO in the future in might be nice to save the filename with the random seed or a hash
    # so that we can just add more bootstraps if needed
    filename <- paste0(output,"/bs", my_seed, ".rds")
    log_info(paste0("Saving bootstrap", i, " to ", filename))
    saveRDS(model, file = filename)

    # Free up memory
    #rm(model)
    
    # Next seed: it's not perfect but just increment the seed
    my_seed <- my_seed +1
    set.seed(my_seed)

    my_filenames[[i]] <- filename
  }

  return(my_filenames)

}

