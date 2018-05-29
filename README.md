# lassoBag
[LassoBag](https://github.com/likelet/lassoBag) provides an implementation of a CART frame for selecting markers from multiple models. A main application of this package is to screen limitted number of variables that are less dependent to
train dataset. Basically, this packages was initially deveploped for adjust LASSO selected results from bootstrapped sample set. Variables with the highest frequency among the several
selected result were considered as stable variables for differ sample set. However, it is usually hard to determine the cutoff in terms of  frequency when applyed in a real dataset. In this package, we introduced 
an permutation test to obtain P values of variables. In addition, it supports several multiple correction methods for  P value adjustment.

#  Installation 
* via github
```R
install.packages("devtools")
devtools::install_github("likselet/lassoBag")
```

