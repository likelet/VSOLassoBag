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
    *the size of resample sample*
    
.. confval:: imputeN

    :type: numeric
    :default: ``1000``
    *the initial permutation times*

.. confval:: imputeN.max

    :type: numeric
    :default: ``2000``
    *the max permutation times. Regardless of whether p has meet the requirement*

.. confval:: permut.increase

    :type: numeric
    :default: ``1000``
    *if the initial imputeN times of permutation doesn't meet the requirement, then we add ¡®permut.increase times of permutation¡¯ to get more random/permutation values*

.. confval:: boot.rep

    :type: bool
    :default: ``TRUE``
    *whether :"sampling with return" or not*

.. confval:: a.family

    :type: string
    :default: ``none``
    *what kind of regression method to use, it should match the type of out.mat*



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
