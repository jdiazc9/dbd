# Code for simulations &amp; analysis in 'Diversity begets diversity under microbial niche construction'
[https://doi.org/10.1101/2022.02.13.480281](https://doi.org/10.1101/2022.02.13.480281)

**./community-simulator** - This directory contains the source files for the community simulator package updated with the functionalities described [here](https://doi.org/10.1073/pnas.2111261119).

**dbd_corr_vs_params.py** - Main simulation script. Runs 100 sets of simulations randomizing model parameters each time, and outputs 3 files:

* dbd_corr_vs_params.txt - Contains information on each set of model parameters used and the diversity slope observed in each case.
* dbd_rpos.txt - Contains the number of non-focal families and focal species in each of the communities for a specific simulation where a positive diversity slope was observed.
* dbd_rneg.txt - Contains the number of non-focal families and focal species in each of the communities for a specific simulation where a negative diversity slope was observed.

**dbd_corr_vs_params.R** - Loads the files above and generates the corresponding plots.

**null_model.R** - An implementation of the null model of community assembly described in the Supplementary Text. Generates the plots in Figure S1.

## Contact

Please direct any questions to [Juan Diaz-Colunga](mailto:juan.diazcolunga@yale.edu).