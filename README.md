# MALVA EXPERIMENTS

This repository contains the scripts to replicate the experiments
described in the MALVA paper.  To create a fully reproducible data
analysis, we used the Snakemake workflow manager.  Moreover, to
simplify the installation of all the required tools, we used conda (we
provide the conda environments).

To clone the repository:
```
git clone --recursive https://github.com/AlgoLab/malva_experiments.git
```

### Prerequisites

- we assume `snakemake` and `conda` to be installed and in your `PATH`

- we used a modified version of `vargeno` (provided as submodule), to
  compile it:
```
cd [malva_experiments_repo]/vargeno
git checkout fix_fasta_header
bash install.sh
```

- while running `discoSNP++`, the output files will be temporarily
  stored in the current directory: we assume there is enough space


### Experiments

Change the `root` folder in the `config.yaml` file with the desired
folder and run:
```
snakemake [-n] --use-conda [-j 4]
```

If you have any problem running the experiments, please contact [Luca
Denti](https://github.com/ldenti).
