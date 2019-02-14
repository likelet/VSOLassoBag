Parameters and Configuration
====


System parameters
----

.. confval:: -profile

    :type: string
    :default: ``none``
    *using predefined profile in your ``conf`` folder*

Mandatery parameters 
----

.. confval:: --reads

    :type: string
    :default: ``none``
    *need to specified to parsing input reads*

.. confval:: --designfile

    :type: string
    :default: ``false``
    *optional when need to perform comparison between conditions* 

.. confval:: --comparefile

    :type: string
    :default: ``false``
    *optional when need to perform comparison between conditions* 


optional parameters 
----

.. confval:: --singleEnd

    :type: boolean
    :default: ``false``
    *set ``true`` to run analysis on single end data set * 

.. confval:: --skip_qc

    :type: boolean
    :default: ``false``
    *set ``true`` to skip fastp processing step * 

.. confval:: --strand

    :type: boolean
    :default: ``false``
    *set ``true`` to run analysis in strand specific mode * 

.. confval:: --skip_multiqc

    :type: boolean
    :default: ``false``
    *set ``true`` to skip multiqc report section from multisamples * 

.. confval:: --without_replicate

    :type: boolean
    :default: ``false``
    *set ``true`` to perform comparison without replicate, using a poisson sourced test instead* 

.. confval:: --skip_gsea

    :type: boolean
    :default: ``false``
    *set ``true`` to skip Gene Set Enrichment Analysis step * 
    
