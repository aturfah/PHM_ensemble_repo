# PHM-Ensemble Analysis Repository

This repo contains code and instructions to replicate the analyses in [Turfah and Wen](https://arxiv.org/html/2611.32945v1). 


## Installing the PHM package

All R code in this repository uses the [PHM package](https://github.com/aturfah/PHM).
Install it from GitHub. It requires R 3.5 or later.

```r
install.packages("remotes")
remotes::install_github("aturfah/PHM")
```

This also installs the CRAN packages PHM depends on (`mclust`, `mvtnorm`,
`ggplot2` and others).

To install from a local clone of the PHM repository instead, run this in a
terminal:

```sh
R CMD INSTALL /path/to/PHM
```

