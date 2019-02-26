Parameters and Configuration
====


Lasso.bag
----

.. confval:: mat

    :type: data.frame
    :default: ``none``
    *independent variables*

.. confval:: out.mat

    :type: data.frame
    :default: ``none``
    *dependent variables, which contains one column or two columns*

.. confval:: bootN

    :type: numeric
    :default: ``1000``
    *the size of resampled sample, only valid when permutation set to TRUE*
    
.. confval:: imputeN

    :type: numeric
    :default: ``1000``
    *the initial permutation times, only valid when permutation set to TRUE*

.. confval:: imputeN.max

    :type: numeric
    :default: ``2000``
    *the max permutation times. Regardless of whether p has meet the requirement,, only valid when permutation set to TRUE*

.. confval:: permut.increase

    :type: numeric
    :default: ``1000``
    *if the initial imputeN times of permutation doesn't meet the requirement, then we add ¡®permut.increase times of permutation¡¯ to get more random/permutation values, only valid when permutation set to TRUE*

.. confval:: boot.rep

    :type: bool
    :default: ``TRUE``
    *whether :"sampling with return" or not, only valid when permutation set to TRUE*

.. confval:: a.family

    :type: string
    :default: ``none``
    *what kind of regression method to use, it should match the type of out.mat*

.. confval:: parallel

    :type: bool
    :default: ``FALSE``
    *whether the script run in parallel, you need to set n.cores in case this package conquers all your cpu resource*

.. confval:: fit.pareto

    :type: string
    :default: ``mle``
    *the method of fitting Generalized Pareto Distribution, alternative choice is "gd", for gradient descend*

.. confval:: permutation

    :type: string
    :default: ``TRUE``
    *to decide whether to do permutation test, if set TRUE, no p value returns*

.. confval:: n.cores

    :type: numeric
    :default: ``1``
    *how many cores/process to be assigned for this function, in Windows, you have to set this to 1*

.. confval:: rd.seed

    :type: numeric
    :default: ``89757``
    *it is the random seed of this function, in case some of the experiments need to be reappeared*

.. confval:: plot.freq

    :type: string/bool
    :default: ``FALSE``
    *whether to show all the non-zero frequency in the final barplot or not. If "full", all the features(including zero frequency) will be plotted. If "part", all the non-zero features will be plotted. If "not", will not print the plot.*

.. confval:: plot.out

    :type: string/bool
    :default: ``FALSE``
    *the path or file's name to save the plot. If set to FALSE, no plot will be output. Default to FALSE*


LessPermutation 
----

.. confval:: X

    :type: vector
    :default: ``none``
    *a union of input data, e.g. c(1,2,3,4,5,6).*

.. confval:: x0

    :type: numeric
    :default: ``none``
    *the observed value* 

.. confval:: fitting.method

    :type: string
    :default: ``gd``
    *the fitting method of General Pareto Distribution(GPD)* 

.. confval:: search.step

    :type: numeric
    :default: ``0.01``
    *the length of step (this param * length(X)) to find threshold, default 0.01* 

.. confval:: fit.cutoff

    :type: numeric
    :default: ``0.05``
    *the cutoff of p value to judge whether it fits well to GPD, default 0.05* 
 
.. confval:: when.to.fit

    :type: numeric
    :default: ``0.05``
    *a cutoff to tell how many sample values are bigger than the observed value then we don't need to fit GPD. it is a portion.* 
