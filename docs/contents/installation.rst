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

You need to input dependent and independent variables and the function will return a data frame which contains the p value and adjusted p value for each features in the dependent variables. The details of its parameters are shown in parameters section.

.. code-block:: bash

    require(glmnet)
    require(POT)
    require(parallel)
    require(ggplot2)
    
    df <- df.test # this is the integrated data of this package
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
    rownames(df) <- df$ID
    df <- df[,which(colnames(df)!="ID")]

    x <- df[,which(!colnames(df) %in% c("Sex","Age","Osstatus","DFSstatus","OS","DFS","LeftOrRight","LymStatus","NI","VI","Stage","remove","riskscore","riskscoreStatus","ageStatus"))]
    
    # cox
    y <- df[,which(colnames(df) %in% c("Osstatus","OS"))]
    y <- y[,c(2,1)]
    colnames(y) <- c("time","status")
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "cox",parallel=F)

    # binomial
    y <- df$riskscoreStatus
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "binomial",parallel=F)

    # gaussian
    y <- df$riskscore
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F)
    
    # if you don't need the result of permutation
    m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F,permutation=FALSE)


2.2 LessPermutation
^^^^^^^^^^^^^^^

You need to input a union of number and an observed number, the function will return the p value of permutation test.

.. code-block:: bash

    x <- rgpd(200, 1, 2, 0.25)
    LessPermutation(x,1,fitting.method='gd')
