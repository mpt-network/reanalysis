
Contains all data and script for the first MPT reanalysis (i.e., which contains between-subjects comparisons only) across several data sets.

**Note:** Many of the binary `R` files (usually saved as `.RData`) require `R` version `3.5.0` or greater.

# `R` and `Rmd` Files 

The following files and folders are contained in this repository:

- `combine_data.R`: COmbines all individual data files into one combined data file `combined_results.RData`
- `fun_prep.R`: Functions for preparing and combining the individual results files. Used by the data preparation scripts.
- `clean_data_files.R`: Resaves all data files in `data` with high compression and removes `data` attributes which could hold the individual-level data such as no data should be part of this repository accidentally. Was run once before comitting data, does not need to be run again.



# Folders

- `data`: Cleaned data from each subgroup. Some of the folders contain other files with further information. 
