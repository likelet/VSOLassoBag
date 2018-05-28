# lassoBag
[LassoBag](https://github.com/likelet/lassoBag) provides an implementation of a CART frame for selecting markers from multiple models. The main application of this package is to screen limited number of variables that are less dependent to
train dataset. Basically, this packages was initially deveploped for adjust LASSO selected results from bootstrap sample sets. Variables with the highest frequency among the several
selected result were considered as stable variables for differ sample set. However, it usually hard to determine the cut frequency when applying in a real dataset. In this package, we intrudoced
an imputation strategy to obtain P values of variables. In addition, it supported sseveral multiple correction method for  P value adjustment.

#  Installation 
* via github
```R
install.packages("devtools")
devtools::install_github("likselet/lassoBag")
```

