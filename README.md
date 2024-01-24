# source-attribution-MLST

This repo describes the source attributions project from
the Enterics Diseases Epidemiology Branch (EDEB)
and the Enteric Diseases Laboratory Branch (EDLB)
in the Division of Foodborne, Waterborne, and Environmental Diseases
at CDC.

The manuscript can be found at <https://doi.org/10.1089/fpd.2023.0046>.

## Installation

Requires R and `devtools`

```shell
devtools::install_github("ncezid-biome/source-attribution-MLST")
```

To install the executable scripts, locate where they were installed,
e.g., `$HOME/R/4/3.2/sourcerer`, and update your `PATH`.
In Linux, it looks like this:

```shell
# Find the installation
find ~ -type d -name sourcerer
# => assume the path is $HOME/R/4/3.2/sourcerer for this example
export PATH=$PATH:$HOME/R/4/3.2/sourcerer/exec
# Double check the installation with the which command
which bootstrapRF.R
which predictRF.R

```


## Usage

You can use this package either with the executable Rscripts
or via the API.

Either way, first, create a model or models with `bootstrapping()` or  `Rscript bootstrapRF.R`.
Then, query the model with `prediction()` or `Rscript predictRF.R`.

Creating a model is done with a "bootstrap" because it represents only one
outcome out of all the stochatic outcomes in the random forest model.
Therefore, the usage is to create many models and then run prediction
on many models to get an aggregated outcome.

### Examples with API

```R
library("sourcerer")

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

# On second thought, get 10 bootstraps instead of the default of 1
opt$bootstraps <- 10

# Slow RF modeling step
rf_filenames <- bootstrapping(opt)


# Start on the prediction step
option_list <- list(
    make_option(c("-m", "--model"), type = "character", help = "A single random forest model RDS file"),
    make_option(c("-q", "--query"), type = "character", help = "A CSV file with two rows: a header and values for an MLST profile. The header should only have columns with relevant loci and not even an identifier for the genome."),
    make_option(c("-t", "--threads"), type = "integer", help = "How many cores to use. Default: 1", default = 1)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Change the model filenames to the ones you already have
opt$model <- rf_filenames

# required options
required_options <- c("model", "query")
for (o in required_options) {
  if (!(o %in% names(opt))) {
    cat("ERROR: Required option", o, "is missing.\n")
    print_help(opt_parser)
    q(status = 1)
  }
}

pred <- prediction(opt)

# Print the first prediction table
write.table(pred[[1]]$predicted,
            file = stdout(),
            sep  = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE )

```

### Examples with executable scripts

#### Basic example

In this example, we make one random forest model under `results/bs23.rds`, and then we query with the provided example MLST profile `data/example_query.csv`.
If more bootstraps were requested than 1, there would be more files under `results/`.
Each random forest number is given a bootstrap number in the file name, with the random seed used.
Seeds are incremented by 1 if more bootstraps are requested.

    Rscript scripts/bootstrapRF.R --input data/isolates_original_plus_new_dec_1_2021.csv.gz -o results --dependent food --core-loci data/cgMLST_loci.csv --starts-with LMO --bootstraps 1 --threads 8 --seed 23
    Rscript scripts/predictRF.R --query data/example_query.csv -m results/bs23.rds --threads 1

#### Parallelized model creation

You can use `--threads` to use multiple cores in the underlying libraries.
However, you can let the operating system parallelize the calls if you have many models to create.

Using this exact method, we call 4 random seeds with `shuf` and parallelize with `xargs` using 4 processors.
The seed is transferred to the program with `--seed`.
We request 10 models using `--bootstraps` and therefore we are getting 4 x 10 = 40 models.
Formally, it is possible that the seeds are close enough to each other that some output models will override each other.
You can try to avoid that by increasing the range in `shuf` or by creating one-off models to fill in the gaps.

    shuf -i 1-9999 -n 4 | \
      xargs -n 1 -P 4 bash -c '
        Rscript scripts/bootstrapRF.R --input data/isolates_original_plus_new_dec_1_2021.csv.gz -o results/LMO0003 --dependent food --core-loci data/cgMLST_loci.csv --starts-with LMO0003 --bootstraps 10 --threads 1 --seed $0
      '

#### Parallelized predicting

This method reads all the rds model files you have created and sends them to `xargs`.
Each model will give a prediction in tab-delimited format.
At the end as an optional step, we suggest using `column -t` which aligns columns in a terminal.

    \ls results/LMO0003/*.rds | \
      xargs -n 1 -P 8 bash -c '
        Rscript scripts/predictRF.R --query data/example_query.csv -m $0 --threads 1
      ' | \
      column -t

### Complete usage

```text
Usage: scripts/bootstrapRF.R [options]


Options:
        -i INPUT, --input=INPUT
                Spreadsheet describing MLST profiles, in csv or csv.gz format.

        -o OUTPUT, --output=OUTPUT
                Directory of bootstrap random forest models to output

        -d DEPENDENT, --dependent=DEPENDENT
                The dependent variable in the spreadsheet. Default:food

        -c CORE-LOCI, --core-loci=CORE-LOCI
                A comma-separated list of core loci to help remove duplicate isolates. These loci must be present as headers in the spreadsheet from --input.

        --starts-with=STARTS-WITH
                The prefix of all independent variables. Default:LMO

        --seed=SEED
                Random seed. Default:23

        -b BOOTSTRAPS, --bootstraps=BOOTSTRAPS
                How many random forest bootstraps to output

        -t THREADS, --threads=THREADS
                How many cores to use. Default:1

        -h, --help
                Show this help message and exit
```

```text
Usage: scripts/predictRF.R [options]


Options:
        -m MODEL, --model=MODEL
                A single random forest model RDS file

        -q QUERY, --query=QUERY
                A CSV file with two rows: a header and values for an MLST profile. The header should only have columns with relevant loci and not even an identifier for the genome.

        -t THREADS, --threads=THREADS
                How many cores to use. Default: 1

        -h, --help
                Show this help message and exit
```

