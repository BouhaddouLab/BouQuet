# BouQuet: Bouhaddou Lab Quantitative proteomics Toolkit

Post-search analysis of Spectronaut DIA data (protein abundance, phosphorylation, ubiquitination):
QC, MSstats quantification, pathway enrichment, kinase activity and PhosphoSitePlus regulatory sites.

## Install (once, and again to update)

```r
install.packages("pak")
pak::pak("BouhaddouLab/BouQuet")
```

This repository is private: you need access to it, and a GitHub token saved in R once
(`usethis::create_github_token()`, then `gitcreds::gitcreds_set()`).

## Start a new analysis

In RStudio: **File > New File > R Markdown > From Template > BouQuet proteomics analysis**.
This makes a new folder with a fresh copy of the notebook. Put your Spectronaut report, `keys.csv`
and `contrasts.txt` in a `data/` folder next to it, edit the Setup chunk, and run the chunks in order.

Without RStudio: `rmarkdown::draft("my_experiment", template = "bouquet", package = "BouQuet")`

## What's in the package

- `R/BouQuet_functions.R`: every function, grouped by step (setup, databases, MSstats helpers, gene IDs, enrichment, kinases)
- `inst/rmarkdown/templates/bouquet/skeleton/skeleton.Rmd`: the notebook
- `inst/extdata/`: built-in tables (PhosphoSitePlus regulatory sites, human kinase list, kinase-substrate list); see the README there

## For maintainers: making a change

1. Edit the function file or the notebook template.
2. Raise `Version:` in `DESCRIPTION` (e.g. 4.0.0 to 4.0.1).
3. Commit and push. Lab members get the change with `pak::pak("BouhaddouLab/BouQuet")`.
