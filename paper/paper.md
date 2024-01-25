---
title: 'Sourcerer: an R package for predicting food vehicles from pathogenic bacterial genomes'
tags:
  - command line
  - random forest
  - source attribution
  - foodborne disease
authors:
  - name: Zhaohui Cui
    affiliation: 1
  - name: Sean D. Browning
    affiliation: 2
    orcid: 0000-0002-4580-1843
  - name: Lee S. Katz
    affiliation: 1
    orcid: 0000-0002-2533-9161
affiliations:
  - name: Division of Foodborne, Waterborne, and Environmental Diseases (DFWED), National Center for Emerging and Zoonotic Infectious Diseases (NCEZID), Centers for Disease Control and Prevention, Atlanta, GA, United States of America
    index: 1
  - name: Office of the Director, National Center for Emerging and Zoonotic Infectious Diseases (NCEZID), Centers for Disease Control and Prevention, Atlanta, GA, United States of America
    index: 2
bibliography: paper.bib
---

## Statement of need

Every year, the United States experiences hundreds or thousands of
foodborne outbreaks [@dewey2016foodborne].
However, a large percentage of bacterial foodborne outbreaks are
caused by unknown etiology [@10.1093/cid/ciab771].
With the advent of real time genome sequencing surveillance,
genomes for foodborne outbreaks are readily available [@jackson2016implementation].
It has been shown that using a random forest (RF) model
can help lend some predictability to finding food sources [@gu2023predicting].
Therefore we have created a RF package in the R language
to help predict food sources from pathogenic bacterial genomes.

## Implementation

We divide the implementation into two parts:
creation of the RF model and
predicting from the RF model with a query.
Both steps have the `randomForestSRC` package at
the core [@ishwaran2022randomforestsrc].

With a spreadsheet of MLST profiles and observed food sources,
a RF model is created with `bootstrapping()`.
Because RF is a stochastic process, it is good practice
to create many RF models with different seeds.

Next, `prediction()` is called with a query MLST profile
and one of the RF models.
This yields percentages of how confident the model is
that the query comes from any one food source.
The current food sources in the example data are:
dairy, fruit, meat, seafood, and vegetable.
However, a user could potentially have different
input data for the model in the `bootstrapping()` step.

If many bootstrap RF models are available, then `prediction()`
can be called on each one, and the predictions can be
aggregated through, e.g., averages and standard deviations.

The main functions are also implemented in command line
scripts `bootstrapRF.R` and `predictionRF.R`.

## Conclusions

## Acknowledgements

## References

