VSOLassoBag installation and usage
==================================

1 Install VSOLassoBag
---------------------

VSOLassoBag can be run in both Windows system and most POSIX systems. The follow codes should be run in R.

 .. code-block:: R   
    
    ## via github
    install.packages("devtools")
    devtools::install_github("likelet/lassoBag")



2 Usage Example
--------------------

2.1 Lasso.bag
^^^^^^^^^^^^^^^  

Lasso.bag is an one-step function that can be easily utilized for selecting important variables from multiple models inherited from R package glmnet. Several methods (Parametric Statistical Test, Curve Elbow Point Detection and Permutation Test) are provided for the cut-off point decision of the importance measure (i.e. observed selection frequency) of variables.

For a start, you need to input dependent variable(s) as a matrix **out.mat**, independent variables also as a matrix **mat**, and **a.family** for the model used by Lasso as determined by the type of dependent variables. The function will return a data frame which contains the important measure for each features in the dependent variables. The details of its parameters are shown in the parameters section.

You can also tune other parameters to better balance the performance and required resource for analysis. Key parameters include *bootN* for bagging times, and *bagFreq.sigMethod* for the cut-off point decision method.


 .. code-block:: R
    
    # load the example data
    data(simulated_example)
    
    # as shown below, 'a.family' is totally determined by the type of 'out.mat', the dependent variable(s)
    
    # binomial
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$label,bootN=100,a.family="binomial",bagFreq.sigMethod="PST")
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$label,bootN=100,a.family="binomial",bagFreq.sigMethod="CEP")
    
    # gaussian
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$y,bootN=100,a.family="gaussian",bagFreq.sigMethod="PST")
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$y,bootN=100,a.family="gaussian",bagFreq.sigMethod="CEP")
    
    # cox
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$surv,bootN=100,a.family="cox",bagFreq.sigMethod="PST")
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$surv,bootN=100,a.family="cox",bagFreq.sigMethod="CEP")
    
    # multinomial
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$multi.label,bootN=100,a.family="multinomial",bagFreq.sigMethod="PST")
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$multi.label,bootN=100,a.family="multinomial",bagFreq.sigMethod="CEP")
    
    # mgaussian
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$multi.y,bootN=100,a.family="mgaussian",bagFreq.sigMethod="PST")
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$multi.y,bootN=100,a.family="mgaussian",bagFreq.sigMethod="CEP")
    
    # poisson
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$pois,bootN=100,a.family="poisson",bagFreq.sigMethod="PST")
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$pois,bootN=100,a.family="poisson",bagFreq.sigMethod="CEP")
    
    # multi-thread processing is supported if run on a multi-thread, forking-supported platform (detailed see R package 'parallel'), which can significantly accelerate the process
    # you can achieve this by flag 'parallel' to TRUE and set 'n.cores' to an integer larger than 1, depending on the available threads
    # multi-thread processing using 2 threads
    res<-Lasso.bag(mat=test.df$x,out.mat=test.df$label,bootN=100,a.family="binomial",bagFreq.sigMethod="PST",parallel=TRUE,n.cores=2)



2.2 Results
^^^^^^^^^^^^^^^

A list with:
1. the result dataframe contains *variable* with selection frequency >=1 and their *Frequency*, the \"P.value\" and the adjusted p value *P.adjust* of each variable (if set *bagFreq.sigMethod* = \"PST\" or \"PERT\"), or the elbow point indicators \"elbow.point\", while elbow point(s) will be marked with \"\*\" (if set *bagFreq.sigMethod* = \"CEP\");
2. other utility results, including permutation results, the regression model built on LASSOBag results.



3 References
--------------------

 - Park H, Imoto S, Miyano S, 2015. \"Recursive Random Lasso (RRLasso) for Identifying Anti-Cancer Drug Targets\", PLoS ONE 10(11): e0141869. https://doi.org/10.1371/journal.pone.0141869 .
 
 - V\. Satopaa, J. Albrecht, D. Irwin and B. Raghavan, 2011. \"Finding a 'Kneedle' in a Haystack: Detecting Knee Points in System Behavior\", 2011 31st International Conference on Distributed Computing Systems Workshops, pp. 166-171. https://doi.org/10.1109/ICDCSW.2011.20 .

 - Simon, Noah, Jerome Friedman, Trevor Hastie, and Robert Tibshirani. 2011. \"Regularization Paths for Cox’s Proportional Hazards Model via Coordinate Descent.\" Journal of Statistical Software, Articles 39 (5): 1–13. https://doi.org/10.18637/jss.v039.i05 .

 - Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010. \"Regularization Paths for Generalized Linear Models via Coordinate Descent.\" Journal of Statistical Software, Articles 33 (1): 1–22. https://doi.org/10.18637/jss.v033.i01 .
