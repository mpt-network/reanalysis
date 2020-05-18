
Contains all data and script for the first MPT reanalysis (i.e., which contains between-subjects comparisons only) across several data sets.

**Note:** Some of the binary `R` files (usually saved as `.RData`) require `R` version `3.5.0` or greater (file names containing `_OLD` should work with older versions of R). 

# `R` Files 

- `combine_data.R`: Combines all individual data files into one combined data file `combined_results.RData`
- `fun_prep.R`: Functions for preparing and combining the individual results files. Used by the data preparation scripts.
- `fun_analysis.R`: Functions for analysis.
- `clean_data_files.R`: Resaves all data files in `data` with high compression and removes `data` attributes which could hold the individual-level data such as no data should be part of this repository accidentally. Was run once before committing data, does not need to be run again.
- `prepare_pairs_data.R`: Combines the individual data files with the covariates and creates the main data set for the meta-analysis containing all pairs of estimates across methods.
- `analyse_pairs.R` and `pairs_plot_code.R`: Older files that do not contain currently used code (e.g., code used for presenting at previous meetings).


# Folders

- `docs`: Contains the following `html` files that are produced by synonymous `RMarkdown` (`.Rmd`) files that can be directly accessed via the links:
    - [`dataset_overview.html`](https://mpt-network.github.io/reanalysis/dataset_overview.html): Provides an overview of all data sets and fitting methods per data set. Also comes with a pairwise-comparison of all core parameters.
    - [`prepared_pairs_data.html`](https://mpt-network.github.io/reanalysis/prepared_pairs_data.html): Overview of data used for meta-analysis and all covariates.
    - [`univariate_relationships_with_abs_dev.html`](https://mpt-network.github.io/reanalysis/univariate_relationships_with_abs_dev.html): Overview of univariate relationships of all covariates with the absolute deviation (i.e., first step of meta-analysis).  
- `data`: Cleaned data from each subgroup. Some of the folders contain other files with further information. 

# Binary Data Files
- `combined_results.RData`: Contains combined results from all subgroups. Main object is `dm`.
- `all_pairs_core.RData`: All parameter pairs (i.e., comparison of group-level estimates across groups) with all covariates. Main data set for meta-analysis.
- `all_pairs_core_full.RData`: All parameter pairs (i.e., comparison of group-level estimates across groups) plus some additional columns and less structure.

Binary data files containing `_OLD` in the file name should be readable by `R` versions prior to `R 3.5.0`.
