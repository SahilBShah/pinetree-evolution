# pinetree-evolution

*Sahil B. Shah, Adam J. Hockenberry, and Claus O. Wilke*

**Reference:**

Pre-print available at: <link>

## Overview and restrictions

This is an evolutionary program that uses pinetree to evolve phage genomes based on user's specifications. These specifications entail target data related to a user-specified gene expression pattern (.tsv file) and the genome configuration of the phage (.yml file). These can be created using the given python files: random_data_generator.py and create_genome_configurations.py. We have provided multiple gene expression and genome configuration files found in the "target" and "genome_configuration" folders, respectively.

**Genome restrictions**

The genome_configuration files have a few restrictions when creating them.

1) Genes cannot overlap in the current version.
2) Intragenic regions muct be at least 30 nucleotides long (the length of a pre-defined ribosome). However, this parameter can be changed by the user by specifying the ribosome size within the "genome_simulator.py" script.

## Steps to run an evolutionary simulation

1) Open the terminal.
2) Using the terminal, navigate to the folder containing the "evolution.py" file. If in the pinetree-evolution directory type: 
```
./src/python
```
3) In the terminal, the command structure to run the evolutionary program is:
```
"python3 evolution.py [target file name with extension (string)] [yaml file name with extension (string)] [run number (int)] [generations (int)] [number of replicates per generation (int)] [display progress bar (bool)]".
```

  **target files are found in the "pinetree-evolution/data/targets/" folder.**
  
  **YAML files are found in the "pinetree-evolution/data/genome_configurations" folder.**
  
4) Press enter to run.

5) Outputs will be saved in the "results" directory based on the target file's name and the simulations replicate number.

## Examples

Example command to run the program with a progress bar displayed:
```
python3 evolution.py paper_data1.tsv testing.yml 1 5000 10 True
```
Without the progress bar:
```
python3 evolution.py paper_data1.tsv testing.yml 1 5000 10
```

