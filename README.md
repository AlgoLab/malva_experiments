# MALVA EXPERIMENTS

This repository contains the scripts to replicate the experiments
described in the MALVA paper.  To create a fully reproducible data
analysis, we used the Snakemake workflow manager.  Moreover, to
simplify the installation of all the required tools, we used conda (we
provide the conda environments).

### Prerequisites

- we assume `snakemake` and `conda` to be installed and in your `PATH`

- we used a modified version of `vargeno`, so you have to download and
  compile it:
```
git clone https://github.com/AlgoLab/malva_experiments.git
cd malva_experiments
git clone https://github.com/ldenti/vargeno.git
git checkout fix_fasta_header
bash install.sh
```

- while running `discoSNP++`, the output files will be temporarily
  stored in the current directory: we assume there is enough space


### Experiments
```
snakemake [-n] --use-conda -s Snakefile.full
snakemake [-n] --use-conda -s Snakefile.half
snakemake [-n] --use-conda -s Snakefile.analysis
```

If you have any problem running the experiments, please contact Luca
Denti.