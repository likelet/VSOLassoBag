#' LassoBag: a lasso based variable selecting CART framework
#'
#' @docType package
#' @name LassoBag
#' @import glmnet
#' @import limma
#' @import POT
#' @import parallel
#' @import ggplot2

#' @param mat sample matrix that each column represent a variable and rows represent sample data points, all the entries in it should be numeric.
#' @param out.mat vector or dataframe with two columns with the same length as the sample size from `mat`
#' @param base.exist indicates whether we have an existed frequency data frame. We set it to be FALSE under most condition. If you don't know what it means, leave it to be FALSE(default)
#' @param a.family a string vector to determine the data type of out.mat
#' @param bootN the size of resampled sample, only valid when permutation set to TRUE
#' @param imputeN the initial permutation times, only valid when permutation set to TRUE
#' @param imputeN.max the max permutation times. Regardless of whether p has meet the requirement,, only valid when permutation set to TRUE
#' @param permut.increase if the initial imputeN times of permutation doesn't meet the requirement, then we add ‘permut.increase times of permutation??? to get more random/permutation values, only valid when permutation set to TRUE
#' @param boot.rep whether :"sampling with return" or not, only valid when permutation set to TRUE
#' @param parallel whether the script run in parallel, you need to set n.cores in case this package conquers all your cpu resource
#' @param fit.pareto the method of fitting Generalized Pareto Distribution, alternative choice is "gd", for gradient descend
#' @param permutation to decide whether to do permutation test, if set FALSE, no p value returns
#' @param n.cores how many cores/process to be assigned for this function, in Windows, you have to set this to 1
#' @param rd.seed it is the random seed of this function, in case some of the experiments need to be reappeared
#' @param plot.freq whether to show all the non-zero frequency in the final barplot or not. If "full", all the features(including zero frequency) will be plotted. If "part", all the non-zero features will be plotted. If "not", will not print the plot.
#' @param plot.out the path or file's name to save the plot. If set to FALSE, no plot will be output. If you run this function in Linux command line, you don't have to set this param for the plot.freq will output your plot to your current working directory with name "Rplot.pdf".Default to FALSE.
#' @param do.plot whether to print a barplot to show the frequency of the output. This package will show you a barplot when it is TRUE. If it is FALSE, the two plot.* params above will be invalid then.
#' @return  a dataframe that contains the frequency, the p value and the adjusted p value of each feature(if you set permutation=T)
#' @examples
#' NULL
# require(glmnet)
# require(POT)
# require(parallel)
# require(ggplot2)
# require(simctest)
# require(survival)
#
# df <- load(df.test.rda) # this is the integrated data of this package
# # change those improper format in df
# to.numeric1 <- as.character(df$riskscoreStatus)
# to.numeric2 <- as.character(df$LeftOrRight)
# to.numeric3 <- as.character(df$Sex)
# to.numeric1[which(to.numeric1=="Low")] <- 0
# to.numeric1[which(to.numeric1=="High")] <- 1
# to.numeric2[which(to.numeric2=="Left")] <- 0
# to.numeric2[which(to.numeric2=="Right")] <- 1
# to.numeric3[which(to.numeric3=="Female")] <- 0
# to.numeric3[which(to.numeric3=="Male")] <- 1
# df$riskscoreStatus <- to.numeric1
# df$LeftOrRight <- to.numeric2
# df$Sex <- to.numeric3
# rownames(df) <- df$ID
# df <- df[,which(colnames(df)!="ID")]
#
# x <- df[,which(!colnames(df) %in% c("Sex","Age","Osstatus","DFSstatus","OS","DFS","LeftOrRight","LymStatus","NI","VI","Stage","remove","riskscore","riskscoreStatus","ageStatus"))]
#
# # cox
# y <- df[,which(colnames(df) %in% c("Osstatus","OS"))]
# y <- y[,c(2,1)]
# colnames(y) <- c("time","status")
# m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "cox",parallel=F)
#
# # binomial
# y <- df$riskscoreStatus
# m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "binomial",parallel=F)
#
# # gaussian
# y <- df$riskscore
# m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F)
#
# # if you don't need the result of permutation
# m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F,permutation=FALSE)



#   Test Package:              'Cmd + Shift + T'

#' @export

 # library(glmnet)
 # library(parallel)
 # library(POT)
 # library(simctest)
 # library(ggplot2)
 # library(survival)
 # source("LessPermutation.R")
 # source("PreFilter.R")
 # source("simpleEstimation.R")
 # source("forwardSelection.R")

Lasso.bag <- function(mat,out.mat, base.exist=F,bootN=1000,imputeN=1000,imputeN.max=2000,permut.increase=100,boot.rep=TRUE,a.family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),additional.info=NULL,
                      bagFreq.sigMethod="simple estimation",local.min.BIC=FALSE,use.gpd=F, fit.pareto="gd",parallel=FALSE,n.cores=1,rd.seed=89757,nfolds=4,lambda.type="lambda.1se",
                      plot.freq="part",plot.out=FALSE, do.plot=FALSE,output_dir=NA,
                      filter_method="auto",inbag.filter=TRUE,filter_thres_method="traditional fdr",filter_thres_fdr=0.05,filter_rank_cutoff=0.05,FilterInitPermutT=100,FilterIncPermutSt=100,FilterMAXPermutT=2000,filter_result_report=TRUE,report_all=TRUE,
                      post.regression=FALSE,post.LASSO=FALSE,pvalue.cutoff=0.05) {
                      
  # bootN is the size of resample sample
  # mat is independent variable
  # out.mat is dependent variable
  # boot.rep is whether :"sampling with return" or not
  # a.family is what kind of regression method to use, it should match the type of out.mat
  # imputeN is the initial permutation times
  # imputeN.max is the max permutation times. Regardless of whether p has meet the requirement
  # permut.increase is: if the initial imputeN times of permutation doesn't meet the requirement, then we add permut.increase times of permutation to get more random values
  # if a.family == "multinomial", then the number of each class of dependent vars should be more than 1
  # additional.info is only for providing additional info to build a cox model, and thus is only for Cox method (a.family == "cox"); a data.frame with same rows as mat/out.mat
  
  if (!is.na(output_dir)){
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    setwd(output_dir)
  }
  
  ## INIT process
  
  ## set random.seed, parallel process
  ## correct input format
  
  set.seed(rd.seed)
  time1<-0
  time2<-0
  time3<-0
  features<-colnames(mat)
  if (!parallel){
    n.cores<-1
  }
  all.num <- bootN * imputeN
  
  if (is.vector(out.mat)){  # force to be matrix
    out.mat<-as.matrix(out.mat,ncol=1)
  }
  
  if (!is.null(additional.info)){
    if (is.vector(additional.info)){
      additional.info<-as.data.frame(additional.info)
    }
    for (i in 1:ncol(additional.info)){
      # For additional.info for building a cox model, all columns should be type of numeric, otherwise it will be forced to be factor first, then transformed to be numeric.
      if (!is.numeric(additional.info[,i])){
        if (!is.factor(additional.info[,i])){
          additional.info[,i]<-as.factor(additional.info[,i])
        }
        additional.info[,i]<-as.numeric(additional.info[,i])
      }
    }
    additional.info<-as.matrix(additional.info)
  }
  
  # simply judge whether dependent vars has the same size of independent vars
  if(nrow(mat)!=nrow(out.mat)){  # TODO: i think that it should be dataframe for that will be more adjustable
    # and it allows us to store different types of data in one mat/df
    warning("incoporate length of matrix and outVarianbles, row No. not match, plz check your input ")
    return("LASSOBag failed due to incorrect input format, row No. of X and Y not match")
  }
  # check column No. of dependent vars
  if (a.family=="cox" & ncol(out.mat)!=2){
    warning("column No. of Y not equals 2 (as 'time' and 'status')")
    return("LASSOBag failed due to incorrect input format, column No. of Y not match for building a Cox model")
  }
  
  # rownames(out.mat) <- c(1:length(rownames(out.mat)))
  # a simple input judgement and transformation
  if (a.family!="cox") {
    colnames(out.mat) <- "out"
  }else{
    colnames(out.mat) <- c("time","status")  #same order as in cv.glmnet function input
  }

  # to keep the same type of each column, default is force them to be "numeric"
  zeromat <- matrix(0,ncol=ncol(mat),nrow=nrow(mat))
  colnames(zeromat) <- colnames(mat)
  rownames(zeromat) <- rownames(mat)
  for (i in 1:ncol(mat)) {
    zeromat[,i] <- as.numeric(mat[,i])
  }
  mat <- zeromat  # it's necessary to substitute in this way for there will be some mistake in later handling if we don't do this
  
  rm(zeromat)
  gc()

  zeromat <- matrix(0,ncol=ncol(out.mat),nrow=nrow(out.mat))
  colnames(zeromat) <- colnames(out.mat)
  rownames(zeromat) <- rownames(out.mat)
  for (i in 1:ncol(out.mat)) {
    zeromat[,i] <- as.numeric(out.mat[,i])
  }
  out.mat <- zeromat  # it's necessary to substitute in this way for there will be some mistake in later handling if we don't do this

  rm(zeromat)
  gc()
  
  ## INIT end
  
  ## out-bag filter
  if (filter_method!="none" & !inbag.filter){
    o<-filters(mat,out.mat,silent=FALSE, additional.info=additional.info,filter_method=filter_method,a.family=a.family,
               filter_thres_method=filter_thres_method,filter_thres_fdr=filter_thres_fdr,filter_rank_cutoff=filter_rank_cutoff,
               parallel=parallel,n.cores=n.cores,FilterInitPermutT=FilterInitPermutT,FilterIncPermutSt=FilterIncPermutSt,FilterMAXPermutT=FilterMAXPermutT,
               filter_result_report=filter_result_report,report_all=report_all,rd.seed=rd.seed)
    mat<-mat[,o$sel]
    time1<-time1+o$time
  }
  
  # lasso bag function in individual iteration
  # boot for one time.
  boot.once<-function(sampleIndex,permutate.index.list){  # TODO: change here to check the input type and take actions to exact types
    
    effectdata.x<-mat[sampleIndex,]
    effectdata.y<-out.mat[sampleIndex,]
    effectdata.y<-as.matrix(effectdata.y,nrow=nrow(out.mat),ncol=ncol(out.mat))
    if (!is.null(additional.info)){
      additional.info.boot<-additional.info[sampleIndex,]
    }else{
      additional.info.boot<-NULL
    }
    
    ## Apply in-bag filter
    if (filter_method!="none" & inbag.filter){
      if (filter_thres_method=="permutation fdr"){
        filter_thres_method<-"rank"
      }
      o<-filters(effectdata.x,effectdata.y,silent=TRUE, additional.info=additional.info.boot,filter_method=filter_method,a.family=a.family,
                 filter_thres_method=filter_thres_method,filter_thres_fdr=filter_thres_fdr,filter_rank_cutoff=filter_rank_cutoff,
                 parallel=parallel,n.cores=n.cores,FilterInitPermutT=FilterInitPermutT,FilterIncPermutSt=FilterIncPermutSt,FilterMAXPermutT=FilterMAXPermutT,
                 filter_result_report=filter_result_report,report_all=report_all,rd.seed=rd.seed)
      effectdata.x<-effectdata.x[,o$sel]
      time1<-time1+o$time
    }
    
    if (!is.null(additional.info)){
      effectdata.x<-cbind(effectdata.x,additional.info.boot)
    }
    
    # apply function for not permutated or permutated lasso
    boot.indiv<-function(PermutateIndex){
      out<-effectdata.y[PermutateIndex,]
      out<-as.matrix(out,nrow=nrow(out.mat),ncol=ncol(out.mat))
      colnames(out)<-colnames(out.mat)
      
      cv.glmmod<-cv.glmnet(x=effectdata.x, y=out, family = a.family,nfolds = nfolds)
      result<-coef(cv.glmmod, s = lambda.type)  #re-check this change
      return(result@Dimnames[[1]][which(as.numeric(result) != 0)])  # return the selected feature names
    }
  
    if (is.null(permutate.index.list)){  #No permutation applied
      permutate.index.list<-list()
      permutate.index.list[[1]]<-c(1:length(sampleIndex))
    }
    
    if (parallel & length(permutate.index.list)>1){
      selecVlist1 <- mclapply(permutate.index.list, boot.indiv,mc.cores = n.cores,mc.preschedule=FALSE,mc.cleanup=TRUE)  # Slower, but can save some memory?
    } else {
      selecVlist1 <- lapply(permutate.index.list, boot.indiv)
    }
    
    cat(paste(Sys.time(), "One boostrap done.", sep = "--"),'\n')
    return(selecVlist1)
  }
  
  
  bagging<-function(permutate.index.list=NULL){
    marker_candidate <- colnames(mat)
    selecVlist1 <- list()
    
    #index new list for bagging
    index.list.bootonce<-list()
    for(i in 1:bootN){
      sampleindex2 <- sample(1:nrow(mat),1*nrow(mat),rep=boot.rep)  # re-sampling, same size
      index.list.bootonce[[i]]<-sampleindex2
    }
    
    if (parallel & is.null(permutate.index.list)){
      selecVlist1 <- mclapply(index.list.bootonce, boot.once, permutate.index.list=permutate.index.list,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)
    }else{
      selecVlist1 <- lapply(index.list.bootonce, boot.once, permutate.index.list=permutate.index.list)
    }
    
    cat(paste(Sys.time(), "Bagging finished ...", sep = "--"),'\n')
    
    if (is.null(permutate.index.list)){
      # non-permutated bagging results
      out.vec<-rep(0,length(marker_candidate))
      names(out.vec)<-marker_candidate
      tablecount1 <- table(unlist(selecVlist1))
      out.vec[intersect(names(tablecount1), names(out.vec))] <- tablecount1[intersect(names(tablecount1), names(out.vec))]
      return(out.vec)  # the output is "what features have been chosen."
    }else{
      # mixed permutation results, need to be seperated first
      pttimes<-length(permutate.index.list)
      bgtimes<-bootN
      out.list<-list()
      for (i in 1:pttimes){
        out.vec<-rep(0,length(marker_candidate))
        names(out.vec)<-marker_candidate
        for (j in 1:bgtimes){
          out.vec.temp<-rep(0,length(marker_candidate))
          names(out.vec.temp)<-marker_candidate
          unl<-unlist(selecVlist1[[j]][i])
          tablecount1 <- table(unl)
          out.vec.temp[intersect(names(tablecount1), names(out.vec.temp))] <- tablecount1[intersect(names(tablecount1), names(out.vec.temp))]
          out.vec<-out.vec+out.vec.temp
        }
        out.list[[i]]<-out.vec  # each "column" store the result of a single permutation test
      }
      return(out.list)
    }
  }
  
  # get observed value
  if (!base.exist) {
    cat(paste(date(), "", sep=" -- start calculating observed frequency "), '\n')
    sts<-Sys.time()
    observed.fre<-bagging()  # observed frequency. The true freq
    time2<-difftime(Sys.time(),sts,units="min")*n.cores
  } else {
    observed.fre <- "Base.exists"
  }
  
  gc()
  
  # preserve only features selected at least 1 time in the bagging process
  sav<-which(observed.fre>0)
  observed.fre<-observed.fre[sav]
  mat<-mat[,sav]
  
  # get a res.df to show first, in case cutpoint-decision method is not set
  res.df<-data.frame(variate=names(observed.fre),Frequency=observed.fre,ori.region=sav)
  names(res.df$Frequency)<-NULL
  rownames(res.df)<-NULL
  
  if (bagFreq.sigMethod=="permutation test" | bagFreq.sigMethod=="PT") {
    # Permutation Tests
    #construct datalist with permutation
    cat(paste(date(), "", sep=" -- permutate index "), '\n')

    get.permutation <- function(N) {
      # N is how many times to do permutations at this function
      # returns out.df
      index.list<-list()
      for (i in 1:N) {
        temp.index <- sample(1:nrow(mat),nrow(mat),rep=F)
        index.list[[i]]<-temp.index
      }  # index.list is of length=imputeN, each index is the random order of the original matrix

      cat(paste(date(), "", sep=" -- permutating "), '\n')
      # do permutation

      # the origin
      sts<-Sys.time()
      permut.list<-bagging(permutate.index.list=index.list)  ## possible multiprocessing operation is set inside the function boot.once
      time3<-difftime(Sys.time(),sts,units="min")*n.cores+time3

      return(permut.list)
    }

    get.plist <- function(permut.list,N) {
      # permut.list is from permutation
      # N is the total permutation times

      out.df<-do.call(cbind.data.frame, permut.list)  # output as a data frame
      # features saved in df will be shown else the entities will be 0
      cat(paste(date(), "", sep=" -- getting pvalue "), '\n')

      #get permutation pvalue
      pvalue.list<-c()
      add.more <- F
      for(i in 1:length(observed.fre)){
        good.fit<-"not_fit"
        temp.list<-out.df[names(observed.fre)[i],]
        #avoid p=0
        # p.value<-(length(temp.list[temp.list>observed.fre[i]])+1)/N  # how much is bigger than observed value
        observed.value <- observed.fre[i]
        # to use GPD to fit right tail of distribution, but sometimes it will throw error because there is no
        # data that meets the requirement for fitting GPD
        if (use.gpd) {
          # show if the GPD
          tryCatch({p.value <- LessPermutation(as.numeric(as.character(temp.list)),observed.value,fitting.method = fit.pareto)
          good.fit <- adtest.gpd(temp.list,observed.fre[i],fitting.method = fit.pareto)},
                   error=function(e){p.value<-(length(temp.list[temp.list>observed.fre[i]])+1)/N;
                   print(p.value);
                   print("no data is bigger than threshold, we will use traditional p-value")})

        } else {
          good.fit <- "fit_good_enough"
        }
        
        if (!exists("p.value")) {
          p.value<-(length(temp.list[temp.list>observed.fre[i]])+1)/N
        }
        
        if (p.value==0) {
          p.value <- (length(temp.list[temp.list>observed.fre[i]])+1)/N
        }
        #avoid p>1
        p.value<-ifelse(p.value>1,1,p.value)
        pvalue.list<-c(pvalue.list,p.value) # for every feature, there will be a p-value, for later adjustment
        # judge whether we need to add more permutation, using good of fit test(ad test for gpd)
        # only if one feature doesn't meet the requirement then we should add more permutations
        # if (!length(temp.list[temp.list>observed.fre[i]])>=10) {
        #   add.more <- T
        # }
        if (good.fit!="fit_good_enough") {
          add.more <- T
        }
      }
      return(list(add.more, out.df, pvalue.list))
    }

    # TODO: to add GPD here to reduce permutation times.
    # we must first do imputeN times permutation
    permut.list <- get.permutation(imputeN)
    if (!base.exist) {
      judgement <- get.plist(permut.list,imputeN)
      add.more <- judgement[[1]]
      out.table <- judgement[[2]]
      total.imputeN <- imputeN
      # TODO: to add more permutation here
      if (add.more) {
        new.permut.list <-c(permut.list)
        while (add.more & total.imputeN<imputeN.max) {
          new.perm <- get.permutation(permut.increase)
          new.permut.list <- c(new.permut.list,new.perm)  # up to then all the permut.list
          new.judgement <- get.plist(new.permut.list,total.imputeN)
          add.more <- new.judgement[[1]]
          total.imputeN <- total.imputeN + permut.increase
          pvalue.list <- new.judgement[[3]]
        }
        total.imputeN <- total.imputeN + permut.increase
      } else {
        pvalue.list <- judgement[[3]]
      }
      #FDR calculation
      cat(paste(date(), "", sep=" -- Pvalue adjusting "), '\n')
      # FDR.list<-p.adjust(pvalue.list,method = "bonferroni")
      FDR.list <- p.adjust(pvalue.list,method = "fdr")
      res.df<-data.frame(variate=names(observed.fre),Frequency=observed.fre,ori.region=sav,P.value=pvalue.list,P.adjust=FDR.list)
      names(res.df$Frequency)<-NULL
      rownames(res.df)<-NULL
      cat(paste(date(), "", sep=" -- Done"), '\n')
    } else {
      out.table<-do.call(cbind.data.frame, permut.list)
      cat(paste(date(), "", sep=" -- Clean up out.df "), '\n')
      res.df <- "Don't care"
      cat(paste(date(), "", sep=" -- Done"), '\n')
    }

  } else {  # do not permutate and output frequency; get P-value w/o permutation tests; use SE or FS to decide cutoff point
    cat(paste(date(), "", sep=" -- Using Non-permutation method to help determine a cutoff point "), '\n')
    
    if (bagFreq.sigMethod=="simple estimation" | bagFreq.sigMethod=="SE"){
      cat("Using Simple Estimation Method...", '\n')
      # By fitting the observed freq to a binomial distribution model to estimate the significant p-value and therefore decide the cutoff
      se<-simpleEstimation(res.df=res.df,bootN=bootN)
      pvalue.list<-se$pvalue
      pin<-se$pi
      FDR.list <- p.adjust(pvalue.list,method = "fdr")
      res.df$P.value<-pvalue.list
      res.df$P.adjust<-FDR.list
      if (do.plot){
        #theoretical<-data.frame(x=rbinom(n=nrow(res.df),size=bootN,prob=pin))
        pdf("ObservedFreqDistribution.pdf")
        print(ggplot(res.df) +
          #geom_histogram(aes(x=x),data=theoretical,inherit.aes=F,fill="#b3b3b3",color="#666666")+ 
          geom_histogram(aes(x=Frequency),fill="#ff9999",color="#cc0000")  + geom_text(aes(x=max(Frequency)*0.8,y=nrow(res.df)*0.1,label=paste0("estimated probability=\n",round(pin,digits=6)))) +
          theme_bw() +
          theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
          xlab("baggingFreq")+ylab("Count"))
        dev.off()
      }
    }
    
    if (bagFreq.sigMethod=="step-wise forward selection" | bagFreq.sigMethod=="forward selection" | bagFreq.sigMethod=="FS"){
      cat("Using Step-wise Forward Selection Method...", '\n')
      res.df<-forwardSelection(mat=mat, out.mat=out.mat, res.df=res.df, a.family=a.family, additional.info=additional.info,local.min.BIC=local.min.BIC)
      if (do.plot){
        x<-round(c(1:nrow(res.df))/nrow(res.df)*100,digits=2)
        y<-res.df$BIC
        ploty<-data.frame(x=x,y=y)
        ploty<-ploty[which(!is.na(ploty$y)),]
        pdf("BIC.pdf")
        print(ggplot(ploty, aes(x=x,y=y)) +
          geom_step() +
          theme_bw() +
          theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
          xlab("baggingFreq Rank (%)")+ylab("BIC"))
        dev.off()
      }
    }
    out.table <- "You said no permutation"
    cat(paste(date(), "", sep=" -- Done"), '\n')
  }
  
  model<-NULL
  if (post.regression){
    # Optional, build a regression model based on the features selected by LASSOBag with custom-defined (or default) cutoff p-value (or minimum BIC, depend on the bagFreq.sigMethod)
    # the features and their coefficient are reported, including the intercept or the coef of additional info
    res.df$coef<-NA
    if (bagFreq.sigMethod=="step-wise forward selection" | bagFreq.sigMethod=="forward selection" | bagFreq.sigMethod=="FS"){
      post.selected<-c(1:which(res.df$BIC==min(res.df$BIC,na.rm=T)))
    }else{
      post.selected<-which(res.df$P.value<pvalue.cutoff)
    }
    effect.mat<-mat[,post.selected]
    res.df.ps<-res.df[post.selected,]
    res.df.NOTps<-res.df[setdiff(c(1:nrow(res.df)),post.selected),]
    if (!is.null(additional.info)){
      effect.mat<-cbind(effect.mat,additional.info)
    }
    
    if (post.LASSO){
      cv.glm<-cv.glmnet(x=effect.mat, y=out.mat, family = a.family,nfolds = nfolds)
      result<-as.numeric(coef(cv.glm, s = lambda.type)) 
      model<-cv.glm
    }else{
      glm<-glmnet(x=effect.mat, y=out.mat, family = a.family,lambda=0)
      result<-as.numeric(coef(glm))
      model<-glm
    }
    names(result)<-NULL
    
    if (a.family=="cox"){
      if (is.null(additional.info)){
        res.df.ps$coef<-result
      }else{
        empty<-as.data.frame(matrix(NA,ncol=ncol(res.df),nrow=ncol(additional.info)))
        colnames(empty)<-colnames(res.df)
        empty$variate<-colnames(additional.info)
        res.df.ps<-rbind(res.df.ps,empty)
        res.df.ps$coef<-result
      }
    }else{
      empty<-as.data.frame(matrix(NA,ncol=ncol(res.df),nrow=1))
      colnames(empty)<-colnames(res.df)
      empty$variate[1]<-"Intercept"
      res.df.ps<-rbind(empty,res.df.ps,stringsAsFactors=F)
      res.df.ps$coef<-result
    }
    res.df<-rbind(res.df.ps,res.df.NOTps,stringsAsFactors=F)
    res.df<-res.df[order(res.df$Frequency,decreasing=T,na.last=T),]
  }
  
  if (do.plot) {
    plot.df <- res.df[which(!is.na(res.df$Frequency)),]
    if (plot.freq=="full") {
      other.variate<-setdiff(features,res.df$variate)
      empty.df<-as.data.frame(matrix(NA,ncol=ncol(plot.df),nrow=length(other.variate)))
      colnames(empty.df)<-colnames(plot.df)
      empty.df$variate<-other.variate
      empty.df$Frequency<-0
      plot.df<-rbind(plot.df,empty.df)
    } else if (plot.freq!=FALSE & plot.freq!="part" & plot.freq!="full"){
      plot.df <- res.df
      print("Actually you need to set plot.freq correctly, here we will plot all features.")
    }
    if (plot.freq!=FALSE) {
      if (plot.out!=F) {  # for saving files
        pdf(plot.out)
        g<-ggplot(plot.df, aes(reorder(variate, -Frequency), Frequency)) +
           geom_bar(stat = "identity") +
           theme_bw() +
           theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
           xlab(label = "Features")
        if (nrow(plot.df)>=30){
          g<-g+theme(axis.text.x=element_blank())
        }
        print(g)
        dev.off()
      }
      # for plot on the screen
      # print(ggplot(plot.df, aes(reorder(variate, -Frequency), Frequency)) +
              # geom_bar(stat = "identity") +
              # theme_bw() +
              # theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
              # xlab(label = NULL))
    }
  }
  time<-time1+time2+time3
  return(list(results=res.df, permutations=out.table,time=time,model=model))
}
