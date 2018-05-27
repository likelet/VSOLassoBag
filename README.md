# lassoBag
[LassoBag](https://github.com/likelet/lassoBag) provided a implementation of a CART frame for selecting markers from multiple models. The main application of this package is to screen limited number of variables that are less dependent to
train dataset. Basically, this packages was initially deveploped for adjust LASSO selected results from bootstrap sample sets. Variables with the highest frequency among the several
selected result were considered as stable variables for dataset. However, it usually hard to determine the cutoff frequency when applied in real dataset. In this package, we intrudoced
an imputation strategy to get P value of each variables. We also performed multiple correction for each P value.

#  Installation 
* via github
```R
install.packages("devtools")
devtools::install_github("likselet/lassoBag")
```

