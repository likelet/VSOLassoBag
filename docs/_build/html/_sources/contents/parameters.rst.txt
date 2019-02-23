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
    *the size of resample sample, only valid when permutation set to TRUE*
    
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
    *whether the script run in parallel, you need to set mc.core in case this package conquers all your cpu resource*

.. confval:: fit.pareto

    :type: string
    :default: ``mle``
    *the method of fitting Generalized Pareto Distribution, alternative choice is "gd", for gradient descend*

.. confval:: permutation

    :type: string
    :default: ``TRUE``
    *to decide whether to do permutation test, if set TRUE, no p value returns*


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
