.. RNAseqPipe documentation master file, created by
   sphinx-quickstart on Wed Feb 13 23:01:53 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RNAseqPipe's documentation!
======================================

`RNAseqPipe <https://github.com/likelet/RNAseqPipe>`_ (RNAseq Pipeline) aims at data exploration of RNAseq data set. It begins with the raw sequencing data and then following a step of quality control. We recommend our users to use the Illumina generated paired-end/singled-end sequences and the sequencing depth should be more than 10M. The length of the reads should be longer than 50bp 

More information can be found in The RNAseqPipe workflow page, or in the project GitHub wiki.
Please be sure that you have all dependencies software or tools preinstalled in your system. Otherwise, we recommended that users employ docker or singularity containers to run the pipe. 

This wiki includes several tutorials, plz following the step by step tutorial to run RNAseqPipe in your local server or cluster. This tutorial explains how to set the parameters in the ``nextflow.config`` file, and describe the files that will be produced in output, while at this page you can know more about How to read the logs. 


Pipeline Steps
----

The pipeline allows you to choose between running either replicate or without replicates.
Choose between workflows by using ``--without_replicate`` or not(default).

    +----------------------------------------------+---------------------+----------------------+
    | **Step**                                     | **With replicates** | **without_replicate**|  
    +----------------------------------------------+---------------------+----------------------+
    |Raw data QC                                   |Fastp                |Fastp                 |
    +----------------------------------------------+---------------------+----------------------+                 
    |Align Reads                                   |STAR                 |STAR                  |
    +----------------------------------------------+---------------------+----------------------+
    |Alignment QC                                  |Qualimap             |Qualimap              |   
    +----------------------------------------------+---------------------+----------------------+
    |Reads counting                                |RSEM                 |RSEM                  |
    +----------------------------------------------+---------------------+----------------------+
    |Matrix collapses                              |DAtools              |DAtools               |
    +----------------------------------------------+---------------------+----------------------+
    |Differential expression                       |DESeq2               |Poisson Test(DAtools) |
    +----------------------------------------------+---------------------+----------------------+
    |Gene Set enrichment Analysis                  |GSEA                 |-                     |
    +----------------------------------------------+---------------------+----------------------+
    |Summary Report                                |MultiQC              |MultiQC               |
    +----------------------------------------------+---------------------+----------------------+
 
.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   contents/introdction
   contents/installation
   contents/userguide
   contents/parameters
   contents/output
   contents/faq
   contents/licencing
   contents/citations

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
