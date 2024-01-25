#' @title Bootstrapping
#'
#' @param input (character) A spreadsheet describing MLST profiles, in csv or csv.gz format.
#' @param output (character) Output folder for bootstrap random forest models
#' @param core_loci (character) A comma-separated list of core loci to help remove duplicate isolates.
#'  These loci must be present as headers in the spreadsheet from `input`.
#' @param dependent (character, default "food") The dependent variable in `input`
#' @param var_prefix (character, Default: "LMO") The prefix of all independent variables.
#' @param seed (numeric, default: 23) Random seed
#' @param bootstraps (numeric, default: 1) How many random forest bootstrap models to output
#' @param ncores (integer, default: 1) How many cores to use. Default:1
#' 
#' @return my_filenames A list of random forest filenames
#' @export
#'
#' @examples TODO
#' See [inst/bootstrapRF.R]
#'
#' @import cluster
#' @import dplyr
#' @importFrom stats cutree hclust quantile rmultinom sd xtabs
#' @importFrom utils read.csv write.table
#' @importFrom magrittr `%>%`
#' @import ggplot2
bootstrapping <- function(input, output, core_loci, dependent = "food", var_prefix = "LMO", seed = 23, bootstraps = 1, ncores = 1L) {

  logger::log_info(paste0("Running with ", ncores," cores and ", bootstraps, " bootstraps"))

  orgOpt <- options()
  options(rf.cores = ncores, mc.cores = ncores)
  on.exit(options(orgOpt))

  logger::log_info(paste0("Loci start with ", var_prefix))

  logger::log_info(paste0("Getting the MLST profiles from ", input))
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
    select("food", starts_with(var_prefix)) %>%
    mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
    mutate(across(starts_with(var_prefix), ~ as.factor(as.character(.x)))) %>%
    mutate(food = as.factor(as.character(food))) %>%
    as.data.frame() # rfsrc() doesn't work with a tibble

  # Bootstrapping
  if(!dir.exists(output)){
    dir.create(output)
  }

  log_info(paste0("Running bootstraps and saving them to ", output, "/*.rds"))

  # Set the inital seed for RF models

  my_filenames <- list()
  for (i in 1:bootstraps){
    #seed <- sample(1:as.integer(.Machine$integer.max))
    log_info(paste0("Modeling rep ", i, " with seed ", seed, "..."))
    model <- rfsrc(food ~ ., train.df.all, importance = T, seed = seed) #

    # Save intermediate results
    # TODO in the future in might be nice to save the filename with the random seed or a hash
    # so that we can just add more bootstraps if needed
    filename <- paste0(output,"/bs", seed, ".rds")
    log_info(paste0("Saving bootstrap", i, " to ", filename))
    saveRDS(model, file = filename)

    # Free up memory
    #rm(model)
    
    # Next seed: it's not perfect but just increment the seed
    seed <- seed +1
    set.seed(seed)

    my_filenames[[i]] <- filename
  }

  return(my_filenames)

}

