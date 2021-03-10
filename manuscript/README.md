# Notes

This directory contains the relevant code and (partial) data to replicate results found in the manuscript. Please see `README.md` in the base directory for an explanation of how to obtain the full (raw) dataset from the zenodo repository.

1. `code` contains all the code including `python` scripts for doing some basic manipulations, `R` markdown files for creating figures, and `bash` scripts that were used to run replicate simulations.
2. `figure_data` contains the *processed* data that can be used to make the figures according to the `R` markdown files.
3. `figure_output` contains the output from running the code on the relevant figure data. These are mostly sub-panels for the larger manuscript figures
4. `manuscript_figures` contains our final main text and supporting figures.
5. Lastly, as noted, this directory should ultimately contain `manuscript_results` (which must be obtained from zenodo). This (large) dataset contains the complete set of raw data that was ultimately parsed and trimmed down to create the processed data found in `figure_data`.

The `bash` scripts are provided as examples and users should note that they are intended for parallel execution of replicate simulations on a cluster with approximately 100 CPUs. Users should take caution before executing any of these on a laptop/desktop computer. These can be customized to run replicates in various batch sizes, but overall runtimes for all simulation patterns and replicates found within this manuscript (over 3000) will take on the order of several weeks to run (assuming batch sizes in the 25 to 50 range). 
