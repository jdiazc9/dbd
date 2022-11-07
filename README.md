# Code for simulations &amp; analysis in 'Diversity begets diversity under microbial niche construction'
[https://doi.org/10.1101/2022.02.13.480281](https://doi.org/10.1101/2022.02.13.480281)

**./community-simulator** - This directory contains the source files for the community simulator package updated with the functionalities described [here](https://doi.org/10.1073/pnas.2111261119).

**./data** - This directory contains the experimental data.

**./plots** - Plots generated with the scripts under the ./scripts directory are saved here.

**./scripts** - Scripts for the analysis of experimental data and for the consumer-resource model simlations:

* dbd_corr_vs_params_expInoc.py - Runs consumer-resource model simulations, saves results in ./simul_runs and ./simul_params directories.
* dbd_analyzeSimuls.R - analyzes simulations results (loaded from the ./simul_runs and ./simul_params directories).
* null_model.R - An implementation of the null model of community assembly described in the Supplementary Text. Generates the plots in Figure S1.

## Contact

Please direct any questions to [Juan Diaz-Colunga](mailto:juan.diazcolunga@yale.edu).