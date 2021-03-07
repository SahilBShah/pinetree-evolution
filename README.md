# pinetree-evolution

*Sahil B. Shah, Alexis M. Hill, Claus O. Wilke, and Adam J. Hockenberry*

# Project overview

This repository contains the code and (most) of the data required to replicate the analyses found in the manuscript "Generating complex patterns of gene expression without regulatory circuits" (see article pre-print at: <https://www.biorxiv.org/content/10.1101/2020.11.25.398248v1>). Briefly, this manuscript describes an evolutionary program that uses `pinetree` (see:<https://github.com/clauswilke/pinetree>) to evolve phage genomes to match a pre-defined user phenotype (a target set of gene expression time-courses). 

To generate all figures in the manuscript, please go to zenodo repository (ADD FINAL ZENODO LINK) to download output from all runs and unzip into a folder called `manuscript/manuscript_results/`. At this point, you should be able to run the R scripts that reside in the `manuscript/code/r/` directory.

## Requirements

Most packages used are part of standard scientific/numeric python distributions. However, as previously noted this repository relies on the user having previously installed `pinetree` (version 0.3.0 or greater).


# Usage notes

The repository is organized into two distinct sections with separate goals / use cases. 

1. Users looking to interrogate and / or replicate the results found in our manuscript can navigate to `manuscript/`, which contains all the code and data needed to recreate the main figures from processed data.
2. Users looking to run their own evolutionary simulations can find the relevant code for this use case within the `src/python/` folder.

## Input files

The main evolutionary simulation program expects a target gene expression file (`.tsv` format) and a starting genome configuration (`.yml` format). These files can be created using the given python files: `gene_expression_generator.py` and `create_genome_configurations.py`. Additionally, we have provided multiple gene expression and genome configuration example/template files that can be found in the `data/targets/` and `data/genome_configurations/` directories, respectively.

The genome configuration files have a few restrictions that should be noted by potential users:

1. Genes cannot overlap in the current version.
2. Intragenic regions must be at least 30 nucleotides long (the length of a pre-defined ribosome). However, this parameter can be changed by the user by specifying the ribosome size within the `genome_simulator.py` script.
3. To fully explore which gene *general* gene expression patterns can be simulated, we need to create every possible combination of gene orders, at the moment. Future releases may allow for mutations which swap the placement of genes on the genome and thus alleviate this restriction. For now, however simulating all possible orderings of 3 genes for instance requires creating 6 separate genome configurations. 

# Running an evolutionary simulation

1. Open the terminal.
2. Using the terminal, navigate to the folder containing the `evolution.py` file. In relative path terms, this can be found: 
```
cd ./src/python
```
3. From the terminal, the command structure to run the evolutionary program is:
  ```
  python3 evolution.py [target file name with extension (string)] [yaml file name with extension (string)] [run number (int)] [generations (int)] [number of replicates per generation (int)] [display progress bar (bool)]
  ```
4. Press enter to run.
5. Outputs will be saved in the `results` directory based on the target file's name and the simulations replicate number.

**Run-time expectations:**
For a three-gene model, running for 5,000 generations, users can expect a run-time on the order of 3-5 hours per replicate. Replicate simulations can be parallelized in bash scripts for users looking to run hundreds or thousands of replicates. This is why we include the `run number (int)` command line option as results will get written with this run number appended to them. Thus, running multiple parallel replicates with different run numbers will ensure that results are not over-written (as would be the case if multiple replicate simulations run simultaneously with the same `run number (int)`). 

# Recreate large-scale simulations

In order to simulate all expression patterns used in the manuscript, run the bash script found in `manuscript/code/bash/evolution_execute.sh`. This should run 50 replicates of each of the expression pattern files (.tsv) used in the study for a total of 3000 simulations.

## Working examples

Example command to run the program with a progress bar displayed (50 generations):
```
python3 evolution.py expression_dynamic1.tsv starting.yml 1 50 10 True
```
Without the progress bar (also 50 generations):
```
python3 evolution.py expression_dynamic1.tsv starting.yml 1 50 10
```
Example command to run large-scale simulations in the background:
```
nohup evolution_execute.sh &
```

# Next steps

Some next steps we have planned to expand the simulator's functionality are:
1. Implementing a swapping genes function so that each gene combination does not have to be simulated. (Currently in development and located in the gene_swaps branch)
2. Convert into a package for ease of access.
