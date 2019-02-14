
A step by step Usage 
====

Step 1 
----
.. note:: Install nextflow and docker as we mentioned above. 

step 2  
----
.. note:: Make config of your own machine or your environment 
 

* Change the ``base.config`` to custmize the avaiable resource of your server/cluster. 

* Modify the ``template.config`` to :   

    * setting the Path of ``DAtools.jar``   (optional if using docker )  

    * setting  the path of RSEM star index and gtf file, a typical command of building index is 

    .. code-block:: bash

        rsem-prepare-reference --gtf your.gtf --star --star-path /user/bin/STAR -p 8 hg19.fa /ref/human_hg19
         
    * setting ths GSEA path if your going to run GSEA analysis in your cases 

    * rename the configure file as you like (e.g. ``myenv.config``)
        
* setting the running option if you want to skip several steps in `nextflow.config` file , and also add the profile under the `profiles section`
   
   * fox example, add the following string after in the profiles

    .. code-block:: bash
    
        myenv {
                includeConfig 'conf/base.config'
                includeConfig 'conf/myenv.config'
            }
        

step 3
----
.. note:: Run analysis with raw data and ``designfile``, ``comparefile`` 



When you have parepared your ``designfile`` and ``comparefile`` properly, then you can start your analysis with the following command  

.. code-block:: bash

    nextflow run pathTo/RNAseqPipe/main.nf -profile env --read "pathToReads/*_{1,2}.fq.gz" --designfile "design.file" --comparefile "compare.txt"
    

then,  just wait for your jod done 

 

