Parameters and Configuration
====


Lasso.bag
----

.. confval:: mat

    :type: matrix
    :default: ``none``
    *independent variables. sample matrix that each column represent a variable and rows represent sample data points, all the entries in it should be numeric.*

.. confval:: out.mat

    :type: vector / data.frame
    :default: ``none``
    *dependent variables, which contains one column or two columns. vector or dataframe with the same rows as the sample size of **mat**.*

.. confval:: observed.fre
    
	:type: data.frame
	:default: ``NULL``
	*data.frame with columns **variable** and **Frequency**, which can be obtained from existed LASSOBag results for re-analysis. A warning will be issued if the variables in **observed.fre** not found in **mat**, and these variables will be excluded.*

.. confval:: bootN
    
	:type: integer
	:default: ``1000``
	*the size of re-sampled samples for bagging; smaller consumes less processing time but may not get robust results.*

.. confval:: boot.rep
    
	:type: boolean
	:default: ``TRUE``
	*whether sampling with return or not in the bagging procedure*

.. confval:: a.family
    
	:type: character
	:default: ``none``
	:options: ``\"gaussian\", \"binomial\", \"poisson\", \"multinomial\", \"cox\", \"mgaussian\"``
	*a character determine the data type of out.mat, the same used in glmnet*

.. confval:: additional.covariable
    
	:type: data.frame
	:default: ``NULL``
	*provide additional covariable(s) to build the cox model, only valid in Cox method (**a.family** == \"cox\"); a data.frame with same rows as **mat***

.. confval:: bagFreq.sigMethod
    
	:type: character
	:default: ``\"CEP\"``
	:options: ``\"CEP\", \"PST\", \"PERT\"``
	*a character to determine the cut-off point decision method for the importance measure (i.e. the observed selection frequency). Supported methods are \"Parametric Statistical Test\" (abbr. \"PST\"), \"Curve Elbow Point Detection\" (\"CEP\") and \"Permutation Test\" (\"PERT\"). The default and preferable method is \"CEP\". The method \"PERT\" is not recommended due to consuming time and memmory requirement*

.. confval:: kneedle.S
    
	:type: numeric
	:default: ``10``
	*numeric, an important parameter that determines how aggressive the elbow points on the curve to be called, smaller means more aggressive and may find more elbow points. Default **kneedle.S**=10 seems fine, but feel free to try other values. The selection of **kneedle.S** should be based on the shape of observed frequency curve. It is suggested to use larger S first*

.. confval:: auto.loose
    
	:type: boolean
	:default: ``TRUE``
	*if TRUE, will reduce **kneedle.S** automatically in case no elbow point is found with the set **kneedle.S**; only valid when **bagFreq.sigMethod** is \"Curve Elbow Point Detection\" (\"CEP\")*

.. confval:: loosing.factor
    
	:type: numeric
	:default: ``0.5``
	*a numeric value range in (0,1), which **kneedle.S** is multiplied by to reduce itself; only valid when **auto.loose** set to TRUE*
	
.. confval:: min.S
    
	:type: numeric
	:default: ``0.1``
	*a numeric value determines the minimal value that **kneedle.S** will be loosed to; only valid when **auto.loose** set to TRUE*
	
.. confval:: use.gpd
    
	:type: boolean
	:default: ``FALSE``
	*whether to fit Generalized Pareto Distribution to the permutation results to accelerate the process. Only valid when **bagFreq.sigMethod** is \"Permutation Test\" (\"PERT\")*
	
.. confval:: fit.pareto
    
	:type: character
	:default: ``gd``
	:options: ``\"gd\", \"mle\"``
	*the method of fitting Generalized Pareto Distribution, default choice is \"gd\", for Gradient Descend, and alternative as \"mle\", for Maximum Likelihood Estimation (only valid in \"PERT\" mode)*

.. confval:: imputeN
    
	:type: integer
	:default: ``1000``
	*the initial permutation times (only valid in \"PERT\" mode)*
	
.. confval:: imputeN.max
    
	:type: integer
	:default: ``2000``
	*the max permutation times. Regardless of whether p-value has meet the requirement (only valid in \"PERT\" mode)*
	
.. confval:: permut.increase
    
	:type: integer
	:default: ``100``
	*if the initial imputeN times of permutation doesn**t meet the requirement, then we add ¡®permut.increase times of permutation to get more random/permutation values (only valid in \"PERT\" mode)*
	
.. confval:: parallel
    
	:type: boolean
	:default: ``FALSE``
	*whether run in parallel mode; you also need to set n.cores to determine how much CPU resource to use*


.. confval:: n.cores
    
	:type: integer
	:default: ``1``
	*how many threads/process to be assigned for this function; more threads used results in more resource of CPU and memory required*

.. confval:: rd.seed
    
	:type: numeric
	:default: ``10867``
	*the random seed of this function, in case some of the experiments need to be reproduced*

.. confval:: nfolds
    
	:type: integer
	:default: ``4``
	*an integer > 2, how many folds to be created for n-folds cross-validation LASSO in cv.glmnet*

.. confval:: lambda.type
    
	:type: character
	:default: ``\"lambda.1se\"``
	:options: ``\"lambda.1se\", \"lambda.min\"``
	*character, which model should be used to obtain the variables selected in one bagging. Default is \"lambda.1se\" for less variables selected and lower probability being over-fitting. See the help of **cv.glmnet** for more details.*

.. confval:: plot.freq
    
	:type: character
	:default: ``\"part\"``
	:options: ``\"part\", \"full\", \"not\"``
	*whether to show all the non-zero frequency in the final barplot or not. If \"full\", all the variables(including zero frequency) will be plotted. If \"part\"(default), all the non-zero variables will be plotted. If \"not\", will not print the plot.*

.. confval:: plot.out
    
	:type: boolean / character
	:default: ``FALSE``
	*the file's name of the frequency plot. If set to FALSE, no plot will be output. If you run this function in Linux command line, you don't have to set this param for the plot.freq will output your plot to your current working directory with name \"Rplot.pdf\".Default to FALSE.*

.. confval:: do.plot
    
	:type: boolean
	:default: ``TRUE``
	*if TRUE generate result plots*

.. confval:: output.dir
    
	:type: character
	:default: ``NA``
	*the path to save result files generated by Lasso.bag (if not existed, will be created). Default is NA, will save in the same space as the current working dir*

.. confval:: filter.method
    
	:type: character
	:default: ``\"auto\"``
	:options: ``\"auto\",\"pearson\", \"spearman\", \"kendall\", \"cox\"``
	*the filter method applied to input matrix; default is \"auto\", automatically select the filter method according to the data type of **out.mat**. Specific supported methods are \"pearson\", \"spearman\", \"kendall\" from **cor.test** method, and \"cox\" from **coxph** method, and \"none\" (no filter applied).*

.. confval:: inbag.filter
    
	:type: boolean
	:default: ``TRUE``
	*if TRUE, apply filters to the re-sampled bagging samples rather than the original samples*

.. confval:: filter.thres.method
    
	:type: character
	:default: ``fdr``
	:options: ``\"fdr\",\"rank\"``
	*the method determines the threshold of importance in filters. Supported methods are \"fdr\" and \"rank\"*

.. confval:: filter.thres.P
    
	:type: numeric
	:default: ``0.05``
	*if **filter.thres.method** is \"fdr\", use **filter.thres.P** as the (adjusted) cut-off p-value*

.. confval:: filter.rank.cutoff
    
	:type: numeric
	:default: ``0.05``
	*if **filter.thres.method** is \"rank\", use **filter.rank.cutoff** as the cut-off rank*

.. confval:: filter.min.variables
    
	:type: integer
	:default: ``-Inf``
	*minimum important variables selected by filters. Useful when building a multi-variable cox model since cox model can only be built on limited variables. Default is -Inf (not applied)*

.. confval:: filter.max.variables
    
	:type: integer
	:default: ``Inf``
	*maximum important variables selected by filters. Useful when building a multi-variable cox model since cox model can only be built on limited variables. Default is Inf (not applied)*

.. confval:: filter.result.report
    
	:type: boolean
	:default: ``TRUE``
	*if TRUE generate filter reports for filter results, only vaild when **inbag.filter** set to FALSE (i.e. only valid in **out-bag** filters mode)*

.. confval:: filter.report.all.variables
    
	:type: boolean
	:default: ``TRUE``
	*if TRUE report all variables in the filter report, only valid when **filter.result.report** set to TRUE*

.. confval:: post.regression
    
	:type: boolean
	:default: ``FALSE``
	*build a regression model based on the variables selected by LASSOBag process*

.. confval:: post.LASSO
    
	:type: boolean
	:default: ``FALSE``
	*build a LASSO regression model based on the variables selected by LASSOBag process, only vaild when **post.regression** set to TRUE*

.. confval:: pvalue.cutoff
    
	:type: numeric
	:default: ``0.05``
	*determine the cut-off p-value for what variables were selected by LASSOBag, only vaild when **post.regression** is TRUE and **bagFreq.sigMethod** set to \"Parametric Statistical Test\" or \"Permutation Test\"*

.. confval:: used.elbow.point
    
	:type: character
	:default: ``\"middle\"``
	:options: ``\"middle\",\"first\",\"last\"``
	*determine which elbow point to use if multiple elbow points were detected for what variables were selected by LASSOBag. Supported methods are \"first\", \"middle\" and \"last\". Default is \"middle\", use the middle one among all detected elbow points. Only vaild when **post.regression** is TRUE and **bagFreq.sigMethod** set to \"Curve Elbow Point Detection\"*

