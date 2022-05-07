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
1. the result dataframe, "results", contains *variable* with selection frequency >=1 and their *Frequency*, the \"P.value\" and the adjusted p value *P.adjust* of each variable (if set *bagFreq.sigMethod* = \"PST\" or \"PERT\"), or the elbow point indicators \"elbow.point\", while elbow point(s) will be marked with \"\*\" (if set *bagFreq.sigMethod* = \"CEP\"). This is the main result VSOLassoBag obtained.
2. other utility results, including permutation results, "permutations", the regression model built on LASSOBag results, "model".

For tutorial purpose, here we used two examples utilizing different cut-off point decision methods to demonstrate how to interpret the returned results.

We used simulated example data for gaussian model from the **simulated_example** from the package for the demonstration.



2.3 Example 1: using "CEP" cut-off point decision method
^^^^^^^^^^^^^^^^

"CEP" (i.e. "Curve Elbow Point Detection") is the default and recommended method for cut-off point decision. Assuming a sharp decreasing of the observed frequency may seperate important features from those unimportant ones, the "CEP" method detects the elbow point(s) on the observed frequency curve, and features with observed frequency higher than the elbow point are inferred important.

There may be more than one elbow point detected on the curve when using loose threshold, so it is recommended to use a stricter threshold first (use a larger *kneedle\.S* ) and auto loose the S parameter in case no elbow point can be found.

The returned result, **res$results**, is a data.frame\:



.. csv-table:: Frozen Delights!
   :header: "Treat", "Quantity", "Description"
   :widths: 15, 10, 30

   "Albatross", 2.99, "On a stick!"
   "Crunchy Frog", 1.49, "If we took the bones out, it wouldn't be
   crunchy, now would it?"
   "Gannet Ripple", 1.99, "On a stick!"



.. csv-table::  
   :widths: "auto"
   :header-rows: 1
   :align: "center"
   
   "variable","Frequency","elbow\.point","Diff","Thres"
   "X_2",100,"",0,0
   "X_7",100,"",0,0
   "X_10",100,"",0,0
   "X_3",99,"",1,0
   "X_6",97,"",2,0
   "X_5",89,"\*",8,3.9426
   "X_9",87,"",2,3.9426
   "X_8",81,"",6,3.9426
   "X_1",60,"\*",21,16.9426
   "X_4",44,"",16,16.9426
   "X_468",27,"\*",17,12.9426
   "X_169",25,"",2,12.9426
   "X_55",19,"\*",6,1.9426
   "X_404",19,"",0,1.9426
   "X_108",18,"",1,1.9426
   "X_265",17,"",1,1.9426
   "X_114",15,"",2,1.9426
   "X_286",15,"",0,1.9426
   "X_236",14,"",1,1.9426
   "X_142",13,"",1,1.9426



(only showing the header and the first 20 rows)

**variable**

The name of the variable.

**Frequency**

The observed frequency of the variable.

**elbow\.point**

Indicator, if detected as an elbow point, it is marked with "\*", otherwise left blank.

**Diff**

The calculated difference.

**Thres**

Threshold, only when the difference is larger than the threshold, it will be detected as an elbow point.

In this example, when using default *kneedle\.S* , 4 elbow points were detected. Generally, one can choose the middle ("median") one as the cut-off point. Here we used the middle one as the cut-off point and obtained 10 important variables (X_2 ~ X_4).

Since X_1 ~ X_10 were set to be important features, the obtained result successfully disrecovered all important features and excluded unimportant ones. 

However, it must be pointed out that in practise, such performance is very **unlikely** to be achieved.







3 References
--------------------

 - Park H, Imoto S, Miyano S, 2015. \"Recursive Random Lasso (RRLasso) for Identifying Anti-Cancer Drug Targets\", PLoS ONE 10(11): e0141869. https://doi.org/10.1371/journal.pone.0141869 .
 
 - V\. Satopaa, J. Albrecht, D. Irwin and B. Raghavan, 2011. \"Finding a 'Kneedle' in a Haystack: Detecting Knee Points in System Behavior\", 2011 31st International Conference on Distributed Computing Systems Workshops, pp. 166-171. https://doi.org/10.1109/ICDCSW.2011.20 .

 - Simon, Noah, Jerome Friedman, Trevor Hastie, and Robert Tibshirani. 2011. \"Regularization Paths for Cox’s Proportional Hazards Model via Coordinate Descent.\" Journal of Statistical Software, Articles 39 (5): 1–13. https://doi.org/10.18637/jss.v039.i05 .

 - Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010. \"Regularization Paths for Generalized Linear Models via Coordinate Descent.\" Journal of Statistical Software, Articles 33 (1): 1–22. https://doi.org/10.18637/jss.v033.i01 .
