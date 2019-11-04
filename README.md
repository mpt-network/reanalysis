
Contains all data and script for the first MPT reanalysis (i.e., which contains between-subjects comparisons only) across several data sets.

**Note:** Many of the binary `R` files (usually saved as `.RData`) require `R` version `3.5.0` or greater.

# `R` Files 

- `combine_data.R`: Combines all individual data files into one combined data file `combined_results.RData`
- `fun_prep.R`: Functions for preparing and combining the individual results files. Used by the data preparation scripts.
- `fun_analysis.R`: Functions for analysis.
- `clean_data_files.R`: Resaves all data files in `data` with high compression and removes `data` attributes which could hold the individual-level data such as no data should be part of this repository accidentally. Was run once before committing data, does not need to be run again.



# Folders

- `docs`: Contains the following `html` files that are produced by synonymous `RMarkdown` (`.Rmd`) files that can be directly accessed via the links:
    - [`dataset_overview.html`](https://mpt-network.github.io/reanalysis/dataset_overview.html): Provides an overview of all data sets and fitting methods per data set. Also comes with a pairwise-comparison of all core parameters.
- `data`: Cleaned data from each subgroup. Some of the folders contain other files with further information. 

# Binary Data Files
- `combined_results.RData`: Contains combined results from all subgroups. Main object is `dm`.
