#' Bootstrap random forest analysis
#' 
#' This function creates one or more random forest models from a spreadsheet of MLST profiles.
#' In many situations, you would want to make 
#' many random forest models and make predictions from
#' each one downstream. 
#' 
#' @param input (character) The input file path for MLST profile data in csv or csv.gz format.
#' @param output (character) The output directory for random forest models.
#' @param ncores (integer, default: 1L) The number of cores to use for parallel processing.
#' @param bootstrap_reps (integer, default: 1L) The number of bootstrap replicates.
#' @param loci_start_with (character, default: "LMO") The prefix for loci names.
#' @param my_seed (integer, default: 23L) The seed for reproducibility.
#' @param ... Additional arguments for future expansion. Currently not used.
#'
#' @return my_filenames A list of random forest filenames
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' # Example usage:
#' model_filenames <- bootstrapping(input = "tests/testthat/isolates_original_plus_new_dec_1_2021.csv.gz", output = "results/",
#'                          ncores = 4, bootstrap_reps = 100)
#' }
#'
#' @import dplyr
#' @importFrom logger log_info
#' @importFrom stats cutree hclust quantile rmultinom sd xtabs
#' @importFrom utils read.csv write.table
#' @importFrom magrittr `%>%`
#' @importFrom randomForestSRC rfsrc
bootstrapping <- function(input, output, ncores = 1L,
                          bootstrap_reps = 1L, loci_start_with = "LMO",
                          my_seed = 23L) {


  if (missing(input)) {
    stop("input is a mandatory argument")
  }
  if (missing(output)){
    stop("output is a mandatory argument")
  }

  log_info(paste0("Running with ",ncores," cores and ", bootstrap_reps, " bootstraps"))

  orgOpt <- options()
  options(rf.cores = ncores, mc.cores = ncores)
  on.exit(options(orgOpt))

  log_info(paste0("Loci start with ", loci_start_with))

  log_info(paste0("Getting the MLST profiles from ", input))
  lm_dat <- read.csv(input, header = TRUE)

  ###################################################################
  #  importance of genes from random forest model based on all genes
  ###################################################################
  log_info(paste0("Filtering input data to remove nearly duplicate profiles and to just view relevant loci"))
  train.df.all <- lm_dat %>%
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
  for (i in seq_len(bootstrap_reps)) {
    #seed <- sample(1:as.integer(.Machine$integer.max))
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

