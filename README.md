
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VSOLassoBag

[VSOLassoBag](https://github.com/likelet/LassoBag) is a variable-selection oriented
LASSO bagging algorithm for **biomarker development in omic-based translational**
**research**, providing a bagging LASSO framework
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

You can install the developed version of VSOLassoBag from github:

``` r
# install.packages("devtools")
devtools::install_github("likelet/VSOLassoBag")
```
`VSOLassoBag` has been released on CRAN, you can also install released version by :

``` r
install.packages("VSOLassoBag")
```

## Documentation

A detailed documentation of VSOLassoBag could be found at [here](https://seqworld.com/VSOLassoBag/)

## Authors

Jiaqi Liang and Chaoye Wang from SYSU

## Supervisor

Qi Zhao from SYSUCC, zhaoqi@sysucc.org.cn

## License 
GNU General Public License V3

## Citation 
Liang J, Wang C, Zhang D, Xie Y, Zeng Y, Li T, Zuo Z, Ren J, Zhao Q. VSOLassoBag: a variable-selection oriented LASSO bagging algorithm for biomarker discovery in omic-based translational research. J Genet Genomics. 2023 Jan 3:S1673-8527(22)00285-5. doi: 10.1016/j.jgg.2022.12.005. Epub ahead of print. PMID: 36608930. [link](https://pubmed.ncbi.nlm.nih.gov/36608930/)
