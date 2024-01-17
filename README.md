# source-attribution-MLST

This repo describes the source attributions project from
the Enterics Diseases Epidemiology Branch (EDEB)
and the Enteric Diseases Laboratory Branch (EDLB)
in the Division of Foodborne, Waterborne, and Environmental Diseases
at CDC.

The manuscript can be found at <https://doi.org/10.1089/fpd.2023.0046>.

## Installation

Requires R version 3.6

There is a conda environment under r-base.yml in this repo to help.

## Usage

First, change the variables in bootstrap/*.R to your liking

* Base directory
* number of cores
* Number of bootstraps
* Maybe other variables at the top of the scripts

Second create the random forest models

    Rscript bootstrapRF.R

Third predict with a profile

    Rscript predictRF.R

