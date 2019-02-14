lassoBag installation and usage
============================

1 Install lassoBag
----------------

lassoBag can be run in both Windows system and most POSIX systems. The follow codes should be run in R.

 .. code-block:: bash   
    
    via github
    install.packages("devtools")
    devtools::install_github("likelet/lassoBag")



2 Usage Example
--------------------

2.1 Lasso.bag
^^^^^^^^^^^^^^^  

You need to input dependent and independent variables and the function will return a data frame which contains the p value and adjusted p value for each features in the dependent variables. The details of its parameters are shown :here:`/parameters`.

.. code-block:: bash

    df <- readRDS("test.rds")
    # change those improper format in df
    to.numeric1 <- as.character(df$riskscoreStatus)
    to.numeric2 <- as.character(df$LeftOrRight)
    to.numeric3 <- as.character(df$Sex)
    to.numeric1[which(to.numeric1=="Low")] <- 0
    to.numeric1[which(to.numeric1=="High")] <- 1
    to.numeric2[which(to.numeric2=="Left")] <- 0
    to.numeric2[which(to.numeric2=="Right")] <- 1
    to.numeric3[which(to.numeric3=="Female")] <- 0
    to.numeric3[which(to.numeric3=="Male")] <- 1
    df$riskscoreStatus <- to.numeric1
    df$LeftOrRight <- to.numeric2
    df$Sex <- to.numeric3

    # cox
    x <- df[,which(!colnames(df) %in% c("Osstatus","OS"))]
    y <- df[,which(colnames(df) %in% c("Osstatus","OS"))]
    y <- y[,c(2,1)]
    colnames(y) <- c("time","status")
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "cox",parallel=F)

    # binomial
    x <- df[,which(colnames(df)!="riskscoreStatus")]
    y <- df$riskscoreStatus
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "binomial",parallel=F)

    # gaussian
    x <- df[,which(colnames(df)!="riskscore")]
    y <- df$riskscore
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F)



To stop nextflow from looking for updates online, you can tell it to run in offline mode by specifying the following environment variable in your ~/.bashrc file:::

.. code-block:: bash

    export NXF_OFFLINE='TRUE'


2.2 Development
^^^^^^^^^^^^^^^

If you would like to make changes to the pipeline, it's best to make a fork on GitHub and then clone the files. Once cloned you can run the pipeline directly as above.


3 Pipeline configuration
------------------------

By default, the pipeline runs with the ``standard`` configuration profile. This uses a number of sensible defaults for process requirements and is suitable for running on a simple (if powerful!) basic server. You can see this configuration in [`conf/base.config`](../conf/base.config).

Be warned of two important points about this default configuration:

1. The default profile uses the ``local`` executor
    * All jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.
    * See the `nextflow docs <https://www.nextflow.io/docs/latest/executor.html`_ for information about running with other hardware backends. Most job scheduler systems are natively supported.
2. Nextflow will expect all software to be installed and available on the `PATH`

3.1 Software deps: Docker
^^^^^^^^^^^^^^^^^^^^^^^^^

First, install docker on your system: `Docker Installation Instructions <https://docs.docker.com/engine/installation/>`_


Then, running the pipeline with the option ``-profile standard,docker`` tells Nextflow to enable Docker for this run. An image containing all of the software requirements will be automatically fetched and used from dockerhub (https://hub.docker.com/r/likelet/RNAseqPipe/).

3.2 Software deps: Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^

If you're not able to use Docker then `Singularity <http://singularity.lbl.gov/>`_ is a great alternative.
The process is very similar: running the pipeline with the option ``-profile standard,singularity`` tells Nextflow to enable singularity for this run. An image containing all of the software requirements will be automatically fetched and used from singularity hub.

If running offline with Singularity, you'll need to download and transfer the Singularity image first:

.. code-block:: bash

    singularity pull --name RNAseqPipe.simg shub://likelet/RNAseqPipe

Once transferred, use ``-with-singularity`` and specify the path to the image file:

.. code-block:: bash

    nextflow run /path/to/circPipe -with-singularity circPipe.simg

Remember to pull updated versions of the singularity image if you update the pipeline.

3.3 Software deps: conda
^^^^^^^^^^^^^^^^^^^^^^^^^

If you're not able to use Docker _or_ Singularity, you can instead use conda to manage the software requirements.
This is slower and less reproducible than the above, but is still better than having to install all requirements yourself!
The pipeline ships with a conda environment file and nextflow has built-in support for this.
To use it first ensure that you have conda installed (we recommend `miniconda<https://conda.io/miniconda.html>`_, then follow the same pattern as above and use the flag ``-profile standard,conda``