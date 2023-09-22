
# Repository for the paper “Survey sampling meets causal inference: a simple method to balance covariates distributions”

# Acknowledgements

This work was financed by the National Science Centre in Poland, OPUS 22
grant no. 2020/39/B/HS4/00941.

# Tutorial for the `jointCalib` package

# Highlights

# Reporduction of the results

## Setup

Install relevant packages for the paper. The `jointCalib` package is
available at CRAN but we install the package from github (development
version). The `IPS` package used in the [Sant’Anna et
al. (2022)](https://onlinelibrary.wiley.com/doi/10.1002/jae.2909) is
available only on github.

``` r
install.packages(c("remotes", "ebal", "mvnfast", "data.table", "ggplot2", "laeken", "xtable", "glue", "stringr"))
remotes::install_github("ncn-foreigners/jointCalib@dev")
remotes::install_github("pedrohcgs/IPS") 
```

## Notebooks and results

- Notebooks:
  - [simulation 1
    (EB)](https://htmlpreview.github.io/?https://raw.githubusercontent.com/ncn-foreigners/paper-note-quantiles-obs-studies/main/codes/1-simulation-eb.html),
  - [simulation 2
    (CBPS)](https://htmlpreview.github.io/?https://raw.githubusercontent.com/ncn-foreigners/paper-note-quantiles-obs-studies/main/codes/1-simulation-ps.html)
  - [processing results for the paper]()
- Simulation results may be found in folder `results/`
