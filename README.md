# G-loop Mining

This repository contains the source code, data, and scripts related to the research on **β-hairpin G-loop motifs** and **minimal G-loop motifs** involved in CRBN neosubstrate recognition. This project predicts and validates CRBN-compatible β-hairpin and helical G-loop proteins, contributing to our understanding of molecular glue degraders and their impact on neosubstrate recruitment. This repository reproduces the results of *Mining the CRBN Target Space Redefines Rules for Molecular Glue-induced Neosubstrate Recognition* (see Reference below).

## Overview

The scripts in this repository enable **mining of G-loop motifs** across protein structures derived from PDB and AlphaFold2 models, based on the β-hairpin G-loop around glycine in CK1α. Identified motifs are further processed and validated for CRBN engagement using various biochemical assays.

### Key Features

- **Parallelized mining**: Efficiently process thousands of structural motifs using the provided Bash script, which parallelizes Python script execution across multiple cores.
- **Motif validation**: Predictions are filtered based on structural compatibility with CRBN.
- **Datasets**: Minimal test datasets are included in this repository. The full list of PDB structures are available at zenodo.org/XXX.
  
![G-loop mining](images/gloop_mining_schematic.png)

## Installation

To run this code, clone the repository and install the required dependencies:

```bash
git clone git@github.com:monterosatx/gloop-mining.git
cd gloop-mining
./build_docker.sh
```

## Running the code

```
./run_docker.sh
cd {/your/gloop/directory}
./run_test.sh
```


## Running the mining for the full proteome. 

First download the AF2 and PDB models from zenodo:

Then, run the full list:

```./run_all_list.sh```

# Reference

_*Mining the CRBN Target Space Redefines Rules for Molecular Glue-induced Neosubstrate Recognition*_

*Authors*: G. Petzold, P. Gainza, S. Annunziato, I. Lamberto, P. Trenh, L. A. McAllister, B. DeMarco, L. Schwander, R. D. Bunker, M. Zlotosch, R. SriRamaratnam, S. Gilberto, G. Langousis, E. J. Donckele, C. Quan, V. Strande, G. M. De Donatis, S. B. Alabi, J. Alers, M. Matysik, C. Staehly, A. Dubois, A. Osmont, M. Garskovas, D. Lyon, L. Wiedmer, V. Oleinikovas, R. Lieberherr, N. T. Rubin, D. T. Lam, N. I. Widlund, A. Ritzen, R. M. Caceres, D. Vigil, J. Tsai, O. Wallace, M. Peluso, A. Sadok, A. Paterson, V. Zarayskiy, B. Fasching, D. Bonenfant, M. Warmuth, J. Castle, S. A. Townson

_bioRxiv 2024.10.07.616933; doi: https://doi.org/10.1101/2024.10.07.616933_
