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



2.2 Parameters
^^^^^^^^^^^^^^^

**mat**

sample matrix that each column represent a variable and rows represent sample data points, all the entries in it should be numeric.

**out.mat**

vector or dataframe with the same rows as the sample size of 'mat'.

**observed.fre**

dataframe with columns 'variable' and 'Frequency', which can be obtained from existed LASSOBag results for re-analysis. A warning will be issued if the variables in 'observed.fre' not found in 'mat', and these variables will be excluded.

**bootN**

the size of re-sampled samples for bagging, default 1000; smaller consumes less processing time but may not get robust results.

**boot.rep**

whether sampling with return or not in the bagging procedure

**a.family**

a character determine the data type of out.mat, the same used in glmnet.

**additional.covariable**

provide additional covariable(s) to build the cox model, only valid in Cox method ('a.family' == "cox"); a data.frame with same rows as 'mat'

**bagFreq.sigMethod**

a character to determine the cut-off point decision method for the importance measure (i.e. the observed selection frequency). Supported methods are "Parametric Statistical Test" (abbr. "PST"), "Curve Elbow Point Detection" ("CEP") and "Permutation Test" ("PERT"). The default and preferable method is "CEP". The method "PERT" is not recommended due to consuming time and memmory requirement.

**kneedle.S**

numeric, an important parameter that determines how aggressive the elbow points on the curve to be called, smaller means more aggressive and may find more elbow points. Default 'kneedle.S'=5 seems fine, but feel free to try other values. The selection of 'kneedle.S' should be based on the shape of observed frequency curve. It is suggested to use larger S first.

**auto.loose**

if TRUE, will reduce 'kneedle.S' in case no elbow point is found with the set 'kneedle.S'; only valid when 'bagFreq.sigMethod' is "Curve Elbow Point Detection" ("CEP").

**loosing.factor**

a numeric value range in (0,1), which 'kneedle.S' is multiplied by to reduce itself; only valid when 'auto.loose' set to TRUE.

**min.S**

a numeric value determines the minimal value that 'kneedle.S' will be loosed to; only valid when 'auto.loose' set to TRUE.

**use.gpd**

whether to fit Generalized Pareto Distribution to the permutation results to accelerate the process. Only valid when 'bagFreq.sigMethod' is "Permutation Test" ("PERT").

**fit.pareto**

the method of fitting Generalized Pareto Distribution, alternative choice is "gd", for gradient descend (only valid in "PERT" mode).

**imputeN**

the initial permutation times (only valid in "PERT" mode).

**imputeN.max**

the max permutation times. Regardless of whether p-value has meet the requirement (only valid in "PERT" mode).

**permut.increase**

if the initial imputeN times of permutation doesn't meet the requirement, then we add ‘permut.increase times of permutation to get more random/permutation values (only valid in "PERT" mode).

**parallel**

whether the script run in parallel mode; you also need to set n.cores to determine how much CPU resource to use.

**n.cores**

how many cores/process to be assigned for this function; more cores used results in more resource of CPU and memory used.

**rd.seed**
the random seed of this function, in case some of the experiments need to be reproduced.

**nfolds**

integer > 2, how many folds to be created for n-folds cross-validation LASSO in cv.glmnet.

**lambda.type**

character, which model should be used to obtain the variables selected in one bagging. Default is "lambda.1se" for less variables selected and lower probability being over-fitting. See the help of 'cv.glmnet' for more details.

**plot.freq**

whether to show all the non-zero frequency in the final barplot or not. If "full", all the variables(including zero frequency) will be plotted. If "part", all the non-zero variables will be plotted. If "not", will not print the plot.

plot.out

the file's name of the frequency plot. If set to FALSE, no plot will be output. If you run this function in Linux command line, you don't have to set this param for the plot.freq will output your plot to your current working directory with name "Rplot.pdf".Default to FALSE.

**do.plot**

if TRUE generate result plots.

**output.dir**

the path to save result files generated by Lasso.bag (if not existed, will be created). Default is NA, will save in the same space as the current working dir.

**filter.method**

the filter method applied to input matrix; default is 'auto', automatically select the filter method according to the data type of 'out.mat'. Specific supported methods are "pearson", "spearman", "kendall" from cor.test method, and "cox" from coxph method, and "none" (no filter applied).

**inbag.filter**

if TRUE, apply filters to the re-sampled bagging samples rather than the original samples; default is TRUE.

**filter.thres.method**

the method determines the threshold of importance in filters. Supported methods are "fdr" and "rank".

**filter.thres.P**

if 'filter.thres.method' is "fdr", use 'filter.thres.P' as the (adjusted) cut-off p-value. Default is 0.05.

**filter.rank.cutoff**

if 'filter.thres.method' is "rank", use 'filter.rank.cutoff' as the cut-off rank. Default is 0.05.

**filter.min.variables**

minimum important variables selected by filters. Useful when building a multi-variable cox model since cox model can only be built on limited variables. Default is -Inf (not applied).

**filter.max.variables**

maximum important variables selected by filters. Useful when building a multi-variable cox model since cox model can only be built on limited variables. Default is Inf (not applied).

**filter.result.report**

if TRUE generate filter reports for filter results, only vaild when 'inbag.filter' set to FALSE (i.e. only valid in out-bag filters mode).

**filter.report.all.variables**

if TRUE report all variables in the filter report, only valid when 'filter.result.report' set to TRUE.

**post.regression**

build a regression model based on the variables selected by LASSOBag process. Default is FALSE.

**post.LASSO**

build a LASSO regression model based on the variables selected by LASSOBag process, only vaild when 'post.regression' set to TRUE.

**pvalue.cutoff**

determine the cut-off p-value for what variables were selected by LASSOBag, only vaild when 'post.regression' is TRUE and 'bagFreq.sigMethod' set to "Parametric Statistical Test" or "Permutation Test".

**used.elbow.point**

determine which elbow point to use if multiple elbow points were detected for what variables were selected by LASSOBag. Supported methods are "first", "middle" and "last". Default is "middle", use the middle one among all detected elbow points. Only vaild when 'post.regression' is TRUE and 'bagFreq.sigMethod' set to "Curve Elbow Point Detection".



2.3 Results
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
