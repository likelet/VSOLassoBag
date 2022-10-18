
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VSOLassoBag

[VSOLassoBag](https://github.com/likelet/LassoBag) is a variable-selection oriented
LASSO bagging algorithm for biomarker development in omic-based translational
research, providing a bagging LASSO framework
for selecting important variables from multiple
models. A main application of this package is to screen limitted number
of variables that are less dependent to train dataset. Basically, this
packages was initially deveploped for adjust LASSO selected results from
bootstrapped sample set. Variables with the highest frequency among the
several selected result were considered as stable variables for differ
sample set. However, it is usually hard to determine the cutoff in terms
of frequency when applyed in a real dataset. In this package, we
introduced several methods, namely (1) curve elbow point
detection, (2) parametrical statistical test and (3) permutation test to
help determine the cut-off point for variables selection.

## Installation

You can install the released version of lassoBag from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LassoBag")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("likelet/VSOLassoBag")
```

## Documentation

A detailed documentation of VSOLassoBag could be found at [here](https://seqworld.com/LassoBag/)

## Authors

Jiaqi Liang and Chaoye Wang from SYSU

## Supervisor

Qi Zhao from SYSUCC, zhaoqi@sysucc.org.cn

## License 
GNU General Public License V3
