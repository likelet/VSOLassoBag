#' One-step main function of VSOLassoBag framework
#'
#' An one-step function that can be easily utilized for selecting important variables from multiple models inherited from R package \emph{glmnet}. Several methods (Parametric Statistical Test, Curve Elbow Point Detection and Permutation Test)  are provided for the cut-off point decision of the importance measure (i.e. observed selection frequency) of variables.
#'
#' @import glmnet
#' @import survival
#' @import ggplot2
#' @import POT
#' @import parallel
#' @import pbapply

#'
#' @param ExpressionData ExpressionData is an object constructed by  SummarizedExperiment. It contains the independent variables matrix and outcome variables matrix.
#' @param outcomevariable Variables which must be the column name of the outcome variables matrix were used to construct models.
#' @param observed.fre dataframe with columns 'variable' and 'Frequency', which can be obtained from existed VSOLassoBag results for re-analysis. A warning will be issued if the variables in `observed.fre` not found in `mat`, and these variables will be excluded.
#' @param bootN the size of re-sampled samples for bagging, default 1000; smaller consumes less processing time but may not get robust results.
#' @param boot.rep whether sampling with return or not in the bagging procedure.
#' @param sample.size The sample size in the bagging space, default is 1 (same sample size as the input sample size).
#' @param a.family a character determine the data type of out.mat, the same used in \code{\link[glmnet]{glmnet}}.
#' @param additional.covariable provide additional covariable(s) to build the cox model, only valid in Cox method (`a.family` == "cox"); a data.frame with same rows as `mat`
#' @param bagFreq.sigMethod a character to determine the cut-off point decision method for the importance measure (i.e. the observed selection frequency). Supported methods are "Parametric Statistical Test" (abbr. "PST"), "Curve Elbow Point Detection" ("CEP") and "Permutation Test" ("PERT"). The default and preferable method is "CEP". The method "PERT" is not recommended due to consuming time and memmory requirement.
#' @param kneedle.S numeric, an important parameter that determines how aggressive the elbow points on the curve to be called, smaller means more aggressive and may find more elbow points. Default `kneedle.S`=10 seems fine, but feel free to try other values. The selection of `kneedle.S` should be based on the shape of observed frequency curve. It is suggested to use larger S first.
#' @param auto.loose if TRUE, will reduce `kneedle.S` in case no elbow point is found with the set `kneedle.S`; only valid when `bagFreq.sigMethod` is "Curve Elbow Point Detection" ("CEP").
#' @param loosing.factor a numeric value range in (0,1), which `kneedle.S` is multiplied by to reduce itself; only valid when `auto.loose` set to TRUE.
#' @param min.S a numeric value determines the minimal value that `kneedle.S` will be loosed to; only valid when `auto.loose` set to TRUE.
#' @param use.gpd whether to fit Generalized Pareto Distribution to the permutation results to accelerate the process. Only valid when `bagFreq.sigMethod` is "Permutation Test" ("PERT").
#' @param fit.pareto the method of fitting Generalized Pareto Distribution, default choice is "gd", for gradient descend, and alternative as "mle", for Maximum Likelihood Estimation (only valid in "PERT" mode).
#' @param imputeN the initial permutation times (only valid in "PERT" mode).
#' @param imputeN.max the max permutation times. Regardless of whether p-value has meet the requirement (only valid in "PERT" mode).
#' @param permut.increase if the initial imputeN times of permutation doesn't meet the requirement, then we add ‘permut.increase times of permutation to get more random/permutation values (only valid in "PERT" mode).
#' @param parallel whether the script run in parallel mode; you also need to set n.cores to determine how much CPU resource to use.
#' @param n.cores how many threads/process to be assigned for this function; more threads used results in more resource of CPU and memory used.
#' @param nfolds integer > 2, how many folds to be created for n-folds cross-validation LASSO in \code{\link[glmnet]{cv.glmnet}}.
#' @param lambda.type character, which model should be used to obtain the variables selected in one bagging. Default is "lambda.1se" for less variables selected and lower probability being over-fitting. See the help of `cv.glmnet` for more details.
#' @param plot.freq whether to show all the non-zero frequency in the final barplot or not. If "full", all the variables(including zero frequency) will be plotted. If "part", all the non-zero variables will be plotted. If "not", will not print the plot.
#' @param plot.out the file's name of the frequency plot. If set to FALSE, no plot will be output. If you run this function in Linux command line, you don't have to set this param for the plot.freq will output your plot to your current working directory with name "Rplot.pdf".Default to FALSE.
#' @param do.plot if TRUE generate result plots.
#' @param output.dir the path to save result files generated by \code{\link[VSOLassoBag]{VSOLassoBag}} (if not existed, will be created). Default is NA, will save in the same space as the current working dir.
#' @param filter.method the filter method applied to input matrix; default is `auto`, automatically select the filter method according to the data type of `out.mat`. Specific supported methods are "pearson", "spearman", "kendall" from \code{\link{cor.test}} method, and "cox" from \code{\link{coxph}} method, and "none" (no filter applied).
#' @param inbag.filter if TRUE, apply filters to the re-sampled bagging samples rather than the original samples; default is TRUE.
#' @param filter.thres.method the method determines the threshold of importance in filters. Supported methods are "fdr" and "rank".
#' @param filter.thres.P if `filter.thres.method` is "fdr", use `filter.thres.P` as the (adjusted) cut-off p-value. Default is 0.05.
#' @param filter.rank.cutoff if `filter.thres.method` is "rank", use `filter.rank.cutoff` as the cut-off rank. Default is 0.05.
#' @param filter.min.variables minimum important variables selected by filters. Useful when building a multi-variable cox model since cox model can only be built on limited variables. Default is -Inf (not applied).
#' @param filter.max.variables maximum important variables selected by filters. Useful when building a multi-variable cox model since cox model can only be built on limited variables. Default is Inf (not applied).
#' @param filter.result.report if TRUE generate filter reports for filter results, only vaild when `inbag.filter` set to FALSE (i.e. only valid in out-bag filters mode).
#' @param filter.report.all.variables if TRUE report all variables in the filter report, only valid when `filter.result.report` set to TRUE.
#' @param post.regression build a regression model based on the variables selected by VSOLassoBag process. Default is FALSE.
#' @param post.LASSO build a LASSO regression model based on the variables selected by VSOLassoBag process, only vaild when `post.regression` set to TRUE.
#' @param pvalue.cutoff determine the cut-off p-value for what variables were selected by VSOLassoBag, only vaild when `post.regression` is TRUE and `bagFreq.sigMethod` set to "Parametric Statistical Test" or "Permutation Test".
#' @param used.elbow.point determine which elbow point to use if multiple elbow points were detected for what variables were selected by VSOLassoBag. Supported methods are "first", "middle" and "last". Default is "middle", use the middle one among all detected elbow points. Only vaild when `post.regression` is TRUE and `bagFreq.sigMethod` set to "Curve Elbow Point Detection".
#'
#' @return  A list with (1) the result dataframe, "results", contains "variable" with selection frequency >=1 and their "Frequency", the "P.value" and the adjusted p value "P.adjust" of each variable (if set `bagFreq.sigMethod` = "PST" or "PERT"), or the elbow point indicators "elbow.point", while elbow point(s) will be marked with "*" (if set `bagFreq.sigMethod` = "CEP"). This is the main result VSOLassoBag obtained. (2) other utility results, including permutation results, "permutations", the regression model built on VSOLassoBag results, "model".
#'
#' @seealso \code{\link[glmnet]{glmnet}} and \code{\link[glmnet]{cv.glmnet}} in R package \emph{glmnet}.
#'
#' @export
#'
#' @examples
#' data("ExpressionData")
#' set.seed(19084)
#'
#' # binomial
#' VSOLassoBag(ExpressionData, "y", bootN=1, a.family="binomial", bagFreq.sigMethod="PST", do.plot = FALSE, plot.freq = "not")
#'
#' \donttest   {
#' # gaussian
#' VSOLassoBag(ExpressionData, "y", bootN=2, a.family="gaussian",
#'             bagFreq.sigMethod="PST", do.plot = FALSE, plot.freq = "not")
#' VSOLassoBag(ExpressionData, "y", bootN=2, a.family="gaussian",
#'             bagFreq.sigMethod="CEP", do.plot = FALSE, plot.freq = "not")
#'
#'
#' # cox
#' VSOLassoBag(ExpressionData, c("time","status"), bootN=2,
#'             a.family="cox", bagFreq.sigMethod="PST", do.plot = FALSE,
#'             plot.freq = "not")
#' VSOLassoBag(ExpressionData, c("time","status"), bootN=2, a.family="cox",
#'             bagFreq.sigMethod="CEP", do.plot = FALSE, plot.freq = "not")
#'
#'
#'
#' # mgaussian
#' VSOLassoBag(ExpressionData, c("multi.label.D_1","multi.label.D_2"), bootN=2,
#'             a.family="mgaussian", bagFreq.sigMethod="PST", do.plot = FALSE,
#'             plot.freq = "not")
#' VSOLassoBag(ExpressionData, c("multi.label.D_1","multi.label.D_2"), bootN=2,
#'             a.family="mgaussian", bagFreq.sigMethod="CEP", do.plot = FALSE,
#'             plot.freq = "not")
#'
#' # poisson
#' VSOLassoBag(ExpressionData, "pois", bootN=10, a.family="poisson",
#'             bagFreq.sigMethod="PST", do.plot = FALSE, plot.freq = "not")
#' VSOLassoBag(ExpressionData, "pois", bootN=2, a.family="poisson",
#'             bagFreq.sigMethod="CEP", do.plot = FALSE, plot.freq = "not")
#'
#' # multi-thread processing is supported if run on a multi-thread,
#' # forking-supported platform (detailed see R package 'parallel'),
#' # which can significantly accelerate the process
#' # you can achieve this by flag 'parallel' to TRUE and set 'n.cores' to an
#' # integer larger than 1, depending on the available threads
#' # multi-thread processing using 2 threads
#' VSOLassoBag(ExpressionData, "y", bootN=1000, a.family="binomial",
#'             bagFreq.sigMethod="PST", do.plot = FALSE, plot.freq = "not",
#'             parallel=TRUE,n.cores=1)
#' }

utils::globalVariables(c("coef", "base.exist", "p.adjust", "pdf",
                         "Frequency", "Diff", "Thres", "dev.off",
                         "median", "reorder", "variable", "k_hyb",
                         "sigma_hyb", "stats"))
VSOLassoBag <- function(ExpressionData, outcomevariable, observed.fre=NULL,
                        bootN=1000,boot.rep=TRUE,sample.size=1,
                        a.family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),additional.covariable=NULL,
                        bagFreq.sigMethod="CEP",
                        kneedle.S=10,auto.loose=TRUE,loosing.factor=0.5,min.S=0.1,
                        use.gpd=FALSE, fit.pareto="gd",imputeN=1000,imputeN.max=2000,permut.increase=100,
                        parallel=FALSE,n.cores=1,
                        nfolds=4,lambda.type="lambda.1se",
                        plot.freq="part",plot.out=FALSE, do.plot=TRUE,output.dir=NA,
                        filter.method="auto",inbag.filter=TRUE,filter.thres.method="fdr",filter.thres.P=0.05,filter.rank.cutoff=0.05,filter.min.variables=-Inf,filter.max.variables=Inf,
                        filter.result.report=TRUE,filter.report.all.variables=TRUE,
                        post.regression=FALSE,post.LASSO=FALSE,pvalue.cutoff=0.05,used.elbow.point="middle"
) {
  mat <- t(assay(ExpressionData))
  out.mat <- colData(ExpressionData)
  out.mat <- as.matrix(out.mat[ ,outcomevariable])

  oldwd <- getwd()
  on.exit(setwd(oldwd))
  if (!is.na(output.dir)){
    if (!dir.exists(output.dir)){
      dir.create(output.dir)
    }
    setwd(output.dir)
  }

  ## INIT process

  ## set random.seed, parallel process
  ## correct input format
  sts<-Sys.time()
  time<-0
  parallel_time<-0
  if (is.null(colnames(mat))){  ## assigned each feature a name as 'X_*' if not provided
    colnames(mat)<-paste0("X_",c(1:ncol(mat)))
  }
  features<-colnames(mat)
  if (!parallel){
    n.cores<-1
  }
  all.num <- bootN * imputeN

  if (!is.matrix(out.mat)){  # force to be matrix
    if (is.null(dim(out.mat))){
      out.mat<-as.matrix(out.mat,ncol=1)
    }else{
      out.mat<-as.matrix(out.mat)
    }
  }

  if (!is.null(additional.covariable)){
    if (is.vector(additional.covariable)){
      additional.covariable<-as.data.frame(additional.covariable)
    }
    for (i in 1:ncol(additional.covariable)){
      # For additional.covariable for building a cox model, all columns should be type of numeric, otherwise it will be forced to be factor first, then transformed to be numeric.
      if (!is.numeric(additional.covariable[,i])){
        if (!is.factor(additional.covariable[,i])){
          additional.covariable[,i]<-as.factor(additional.covariable[,i])
        }
        additional.covariable[,i]<-as.numeric(additional.covariable[,i])
      }
    }
    additional.covariable<-as.matrix(additional.covariable)
  }

  # simply judge whether dependent vars has the same size of independent vars
  if(nrow(mat)!=nrow(out.mat)){  # TODO: i think that it should be dataframe for that will be more adjustable
    # and it allows us to store different types of data in one mat/df
    warning("incoporate length of independent variables (X) and dependent variables (Y), samples not match, plz check your input ")
    stop("VSOLassoBag failed due to incorrect input dimension, samples of X and Y not match")
  }
  # check column No. of dependent vars
  if (a.family=="cox" & ncol(out.mat)!=2){
    warning("No. of features of Y not equals 2 (as 'time' and 'status')")
    stop("VSOLassoBag failed due to incorrect input format, No. of features of Y not correct for building a Cox model")
  }

  # a simple input judgement and transformation
  if (a.family!="cox") {    ## assigned each dependent variables a name as 'D_*' if not provided
    if (is.null(colnames(out.mat))){
      colnames(out.mat) <- paste0("D_",c(1:ncol(out.mat)))
    }
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
    if (typeof(out.mat[,i]) %in% c('numeric','factor','integer')){
      zeromat[,i] <- as.numeric(factor(out.mat[,i]))
    }else{
      zeromat[,i] <- as.numeric(out.mat[,i])
    }

  }
  out.mat <- zeromat  # it's necessary to substitute in this way for there will be some mistake in later handling if we don't do this

  rm(zeromat)
  gc()

  if(sum(is.na(out.mat)) > 0){
    warning("No characters allowed in the response variable, please convert to numeric format")
  }

  ## INIT end
  message(paste(date(), "", sep=" -- INIT process completed, start analyzing "), '\n')

  ## out-bag filter
  if (filter.method!="none" & !inbag.filter){
    o<-filters(mat,out.mat, additional.covariable=additional.covariable,filter.method=filter.method,
               filter.thres.method=filter.thres.method,a.family=a.family,filter.thres.P=filter.thres.P,filter.rank.cutoff=filter.rank.cutoff,
               filter.min.variables=filter.min.variables,filter.max.variables=filter.max.variables,
               parallel=parallel,n.cores=n.cores,
               silent=FALSE,filter.result.report=filter.result.report,filter.report.all.variables=filter.report.all.variables)
    mat<-mat[,o$sel]
    message(paste(date(), "", sep=" -- out-bag filter completed "), '\n')
  }

  # lasso bag function in individual iteration
  # boot for one time.
  boot.once<-function(sampleIndex,permutate.index.list,paralleled=FALSE){  # TODO: change here to check the input type and take actions to exact types

    effectdata.x<-mat[sampleIndex,]
    effectdata.y<-out.mat[sampleIndex,]
    effectdata.y<-as.matrix(effectdata.y,nrow=nrow(out.mat),ncol=ncol(out.mat))
    if (!is.null(additional.covariable)){
      additional.covariable.boot<-additional.covariable[sampleIndex,]
    }else{
      additional.covariable.boot<-NULL
    }

    if (parallel){
      if (paralleled){  ## already paralled in upper workerspace, disable parallel here
        sub_parallel<-FALSE
      }else{
        sub_parallel<-TRUE
      }
    }else{
      sub_parallel<-FALSE
    }

    ## Apply in-bag filter
    if (filter.method!="none" & inbag.filter){
      if (filter.thres.method=="permutation fdr"){
        filter.thres.method<-"rank"
      }
      o<-filters(effectdata.x,effectdata.y, additional.covariable=additional.covariable.boot,
                 filter.method=filter.method,a.family=a.family,filter.thres.method=filter.thres.method,filter.thres.P=filter.thres.P,filter.rank.cutoff=filter.rank.cutoff,
                 filter.min.variables=filter.min.variables,filter.max.variables=filter.max.variables,
                 parallel=sub_parallel,n.cores=n.cores,
                 silent=TRUE,filter.result.report=filter.result.report,filter.report.all.variables=filter.report.all.variables)
      effectdata.x<-effectdata.x[,o$sel]
    }
    ##

    if (!is.null(additional.covariable)){
      effectdata.x<-cbind(effectdata.x,additional.covariable.boot)
    }

    # apply function for not permutated or permutated lasso
    boot.indiv<-function(PermutateIndex){
      out<-effectdata.y[PermutateIndex,]
      out<-as.matrix(out,nrow=nrow(out.mat),ncol=ncol(out.mat))
      colnames(out)<-colnames(out.mat)

      cv.glmmod<-try(cv.glmnet(x=effectdata.x, y=out, family = a.family,nfolds = nfolds),silent=TRUE)
      if (is(cv.glmmod, 'try-error')){
        return(character(0))
      }
      result<-coef(cv.glmmod, s = lambda.type)  #re-check this change
      if (a.family=="multinomial" | a.family=="mgaussian"){
        sel<-character(0)
        for (i in c(1:length(result))){
          sel<-union(sel,result[[i]]@Dimnames[[1]][which(as.numeric(result[[i]]) != 0)])
        }
        return(sel)
      }else{
        return(result@Dimnames[[1]][which(as.numeric(result) != 0)])  # return the selected feature names
      }
    }

    if (is.null(permutate.index.list)){  #No permutation applied
      permutate.index.list<-list()
      permutate.index.list[[1]]<-c(1:length(sampleIndex))
    }

    if (parallel & length(permutate.index.list)>1){
      gc()
      parallel_sts<-Sys.time()
      selecVlist1 <- mclapply(permutate.index.list, boot.indiv,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)  # Slower, but can save some memory
      parallel_time<<-parallel_time+difftime(Sys.time(),parallel_sts,units="min")
    } else {
      selecVlist1 <- lapply(permutate.index.list, boot.indiv)
    }
    return(selecVlist1)
  }

  bagging<-function(permutate.index.list=NULL){
    marker_candidate <- colnames(mat)
    selecVlist1 <- list()

    #index new list for bagging
    index.list.bootonce<-list()
    for(i in 1:bootN){
      sampleindex2 <- sample(1:nrow(mat), sample.size*nrow(mat), replace = boot.rep)  # re-sampling, same size
      index.list.bootonce[[i]]<-sampleindex2
    }

    if (parallel & is.null(permutate.index.list)){
      gc()
      parallel_sts<-Sys.time()
      selecVlist1 <- mclapply(index.list.bootonce, boot.once, permutate.index.list=permutate.index.list,paralleled=TRUE,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)
      parallel_time<<-parallel_time+difftime(Sys.time(),parallel_sts,units="min")
    }else{
      selecVlist1 <- pblapply(index.list.bootonce, boot.once, permutate.index.list=permutate.index.list)
    }

    message(paste(Sys.time(), "Bagging finished ...", sep = "--"),'\n')

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

  # get observed frequency
  if (is.null(observed.fre)) {
    message(paste(date(), "", sep=" -- start calculating observed frequency "), '\n')
    observed.fre<-bagging()  # observed frequency. The true freq

    # preserve only features selected at least 1 time in the bagging process
    sav<-which(observed.fre>0)
    observed.fre<-observed.fre[sav]
    #mat<-mat[,sav]

    # get a res.df to show first, in case cutpoint-decision method is not set
    res.df<-data.frame(variable=names(observed.fre),Frequency=observed.fre,ori.region=sav)
    names(res.df$Frequency)<-NULL
    rownames(res.df)<-NULL
  } else {
    message(paste(date(), "", sep=" -- input existed observed frequency "), '\n')
    if (is(observed.fre, 'data.frame') | length(which(!c("variable","Frequency") %in% colnames(observed.fre)))>0){
      warning("data of observed frequency not correct, must be a data.frame with columns 'variable','Frequency', please check your input data", '\n')
      stop("Input observed frequency Error")
    }
    ## read in standard feature_freq.txt output
    ## must be a data.frame with columns "variable","Frequency","ori.region"
    res.df<-observed.fre[,c("variable","Frequency")]
    ord<-as.integer(factor(res.df$variable,levels=features))
    res.df$ori.region<-ord
    if (anyNA(ord)){
      warning("Some variables in observed.fre not found in mat, and these variables will be excluded from the analysis.")
      res.df<-res.df[which(!is.na(res.df$ori.region)),]
    }
  }
  gc()

  if (!(bagFreq.sigMethod %in% c("Permutation Test","PERT","Parametric Statistical Test","PST","Curve Elbow Point Detection","CEP"))){
    if (bagFreq.sigMethod!="none"){
      warning("bagFreq.sigMethod not correctly set, got observed frequency and exit.",'\n')
    }
    res.df$ori.region<-NULL
    return(list(results=res.df, permutations=NULL,model=NULL))
  }

  if (bagFreq.sigMethod=="Permutation Test" | bagFreq.sigMethod=="PERT") {
    # Permutation Tests
    #construct datalist with permutation
    message(paste(date(), "", sep=" -- permutate index "), '\n')

    get.permutation <- function(N) {
      # N is how many times to do permutations at this function
      # returns out.df
      index.list<-list()
      for (i in 1:N) {
        temp.index <- sample(1:nrow(mat), nrow(mat), replace = FALSE)
        index.list[[i]]<-temp.index
      }  # index.list is of length=imputeN, each index is the random order of the original matrix

      message(paste(date(), "", sep=" -- permutating "), '\n')
      # do permutation

      # the origin
      permut.list<-bagging(permutate.index.list=index.list)  ## possible multiprocessing operation is set inside the function boot.once

      return(permut.list)
    }

    get.plist <- function(permut.list,N) {
      # permut.list is from permutation
      # N is the total permutation times

      out.df<-do.call(cbind.data.frame, permut.list)  # output as a data frame
      # features saved in df will be shown else the entities will be 0
      message(paste(date(), "", sep=" -- getting pvalue "), '\n')

      #get permutation pvalue
      pvalue.list<-c()
      add.more <- FALSE
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
          message(p.value);
          message("no data is bigger than threshold, we will use traditional p-value")})
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
        #   add.more <- TRUE
        # }
        if (good.fit!="fit_good_enough") {
          add.more <- TRUE
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
      message(paste(date(), "", sep=" -- Pvalue adjusting "), '\n')
      # FDR.list<-p.adjust(pvalue.list,method = "bonferroni")
      FDR.list <- p.adjust(pvalue.list,method = "fdr")
      res.df<-data.frame(variable=names(observed.fre),Frequency=observed.fre,ori.region=sav,P.value=pvalue.list,P.adjust=FDR.list)
      names(res.df$Frequency)<-NULL
      rownames(res.df)<-NULL
      res.df<-res.df[order(res.df$Frequency,decreasing=TRUE),]
      message(paste(date(), "", sep=" -- Done"), '\n')
    } else {
      out.table<-do.call(cbind.data.frame, permut.list)
      message(paste(date(), "", sep=" -- Clean up out.df "), '\n')
      res.df <- "Don't care"
      message(paste(date(), "", sep=" -- Done"), '\n')
    }

  } else {  # do not permutate and output frequency; get P-value w/o permutation tests; use PST or CEP to decide cutoff point
    message(paste(date(), "", sep=" -- Using Non-permutation method to determine a cutoff point "), '\n')

    if (bagFreq.sigMethod=="Curve Elbow Point Detection" | bagFreq.sigMethod=="CEP"){
      message("Detecting Elbow Points on the Observed Frequency Curve...", '\n')
      # Detect elbow points on the observed frequency curve and consider them as cutoff point
      res.df<-kneedle(res=res.df,S=kneedle.S,auto.loose=auto.loose,loosing.factor=loosing.factor,min.S=min.S)
      if (length(which(res.df$elbow.point=="*"))==0){
        warning("No elbow point detected on the curve, use simple estimation method instead...",'\n')
        bagFreq.sigMethod<-"PST"
      }else{
        if (length(which(res.df$elbow.point=="*"))>5){
          message("INFO: >5 elbow points detected on curve, may need to use larger kneedle.S to get a more conservative result; also refer to the 'ObservedFreqCurve' plot for kneedle.S adjustment.",'\n')
        }
        if (do.plot){
          pdf("ObservedFreqCurve.pdf")
          x<-c(1:nrow(res.df))
          markPoints<-which(res.df$elbow.point=="*")
          print(ggplot(res.df) +
                  geom_step(aes(x=x,y=Frequency),color="black",size=0.5)+
                  geom_line(aes(x=x,y=Diff),color="red",size=0.5)+
                  geom_step(aes(x=x,y=Thres),color="#6666ff",linetype="twodash",size=0.4)+
                  geom_vline(xintercept=markPoints,color="#0000cc",linetype="dashed",size=0.3) +
                  theme_bw() +
                  theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
                  xlab("Variables Frequency Rank\n(black: frequency; dashed vertical: elbow point;\nred: difference; dashed blue: threshold)")+
                  ylab("Frequency"))
          dev.off()
        }
      }
    }

    if (bagFreq.sigMethod=="Parametric Statistical Test" | bagFreq.sigMethod=="PST"){
      message("Using Parametric Statistical Test...", '\n')
      # By fitting the observed freq to a binomial distribution model to estimate the significant p-value for each feature being "important" and therefore decide the cutoff
      se<-simpleEstimation(res.df=res.df,bootN=bootN)
      pvalue.list<-se$pvalue
      pin<-se$pi
      FDR.list <- p.adjust(pvalue.list,method = "fdr")
      res.df$P.value<-pvalue.list
      res.df$P.adjust<-FDR.list
      if (do.plot){
        pdf("ObservedFreqDistribution.pdf")
        print(ggplot(res.df) +
                geom_histogram(aes(x=Frequency),fill="#ff9999",color="#cc0000")  + geom_text(aes(x=max(Frequency)*0.8,y=nrow(res.df)*0.1,label=paste0("average selection ratio=\n",round(pin,digits=4)))) +
                theme_bw() +
                theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
                xlab("Observed Selection Frequency")+ylab("Variables Count"))
        dev.off()
      }
      res.df<-res.df[order(res.df$Frequency,decreasing=TRUE),]
    }

    out.table <- NULL
    message(paste(date(), "", sep=" -- Done"), '\n')
  }

  model<-NULL
  if (post.regression){
    # Optional, build a regression model based on the features selected by VSOLassoBag with custom-defined (or default) cutoff p-value
    # the features and their coefficient are reported, including the intercept or the coef of additional covariable
    message("Building a regression model using features selected by VSOLassoBag...", '\n')
    if (bagFreq.sigMethod=="Permutation Test" | bagFreq.sigMethod=="PERT" | bagFreq.sigMethod=="Parametric Statistical Test" | bagFreq.sigMethod=="PST"){
      post.selected<-res.df$ori.region[which(res.df$P.adjust<pvalue.cutoff)]
    }else{
      if (bagFreq.sigMethod=="Curve Elbow Point Detection" | bagFreq.sigMethod=="CEP"){
        candidates<-which(res.df$elbow.point=="*")
        if (used.elbow.point=="first"){
          post.selected<-res.df$ori.region[c(1:(candidates[1]-1))]
        }else{
          if (used.elbow.point=="last"){
            post.selected<-res.df$ori.region[c(1:(candidates[length(candidates)]-1))]
          }else{
            post.selected<-res.df$ori.region[c(1:(candidates[floor(median(length(candidates)))]-1))]
            if (used.elbow.point!="middle"){
              warning("used.elbow.point not correctly set, will use default method 'middle'.")
            }
          }
        }
      }
    }
    effect.mat<-mat[,post.selected]
    if (!is.null(additional.covariable)){
      effect.mat<-cbind(effect.mat,additional.covariable)
    }
    if (post.LASSO){
      message("Using LASSO model...", '\n')
      glm<-try(cv.glmnet(x=effect.mat, y=out.mat, family = a.family,nfolds = nfolds),silent=TRUE)
    }else{
      glm<-try(glmnet(x=effect.mat, y=out.mat, family = a.family,lambda=0),silent=TRUE)  ## a no-penalty regression
    }
    if (is(glm, 'try-error')){
      warning("Failed to build a regression model, return NULL.")
    }else{
      model<-glm
    }
    message(paste(date(), "", sep=" -- Done"), '\n')
  }

  if (do.plot) {
    plot.df <- res.df[which(!is.na(res.df$Frequency)),]
    if (plot.freq=="full") {
      other.variable<-setdiff(features,res.df$variable)
      empty.df<-as.data.frame(matrix(NA,ncol=ncol(plot.df),nrow=length(other.variable)))
      colnames(empty.df)<-colnames(plot.df)
      empty.df$variable<-other.variable
      empty.df$Frequency<-0
      plot.df<-rbind(plot.df,empty.df)
    } else if (plot.freq!=FALSE & plot.freq!="part" & plot.freq!="full"){
      plot.df <- res.df
      message("Actually you need to set plot.freq correctly, here we will plot all variables.")
    }
    if (plot.freq!=FALSE) {
      if (plot.out!=FALSE) {  # for saving files
        pdf(plot.out)
        gg<-ggplot(plot.df, aes(reorder(variable, -Frequency), Frequency)) +
          geom_bar(stat = "identity") +
          theme_bw() +
          theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
          xlab(label = "Variables")
        if (nrow(plot.df)>=30){
          gg<-gg+theme(axis.text.x=element_blank())
        }
        print(gg)
        dev.off()
      }
    }
  }
  time<-difftime(Sys.time(),sts,units="min")
  res.df$ori.region<-NULL
  return(list(results=res.df, permutations=out.table,model=model))
}
