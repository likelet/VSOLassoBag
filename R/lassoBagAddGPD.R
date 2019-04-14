#' @title LASSO-bagging
#' a lasso based variable selecting CART framework

#' @param mat sample matrix that each column represent a variable and rows represent sample data points, all the entries in it should be numeric.
#' @param out.mat vector or dataframe with two columns with the same length as the sample size from `mat`
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
#' @return  a dataframe that contains the frequency, the p value and the adjusted p value of each feature(if you set permutation=T)
#' @examples
#' require(glmnet)
#' require(POT)
#' require(parallel)
#' require(ggplot2)
#'
#' df <- df.test # this is the integrated data of this package
#' # change those improper format in df
#' to.numeric1 <- as.character(df$riskscoreStatus)
#' to.numeric2 <- as.character(df$LeftOrRight)
#' to.numeric3 <- as.character(df$Sex)
#' to.numeric1[which(to.numeric1=="Low")] <- 0
#' to.numeric1[which(to.numeric1=="High")] <- 1
#' to.numeric2[which(to.numeric2=="Left")] <- 0
#' to.numeric2[which(to.numeric2=="Right")] <- 1
#' to.numeric3[which(to.numeric3=="Female")] <- 0
#' to.numeric3[which(to.numeric3=="Male")] <- 1
#' df$riskscoreStatus <- to.numeric1
#' df$LeftOrRight <- to.numeric2
#' df$Sex <- to.numeric3
#' rownames(df) <- df$ID
#' df <- df[,which(colnames(df)!="ID")]
#'
#' x <- df[,which(!colnames(df) %in% c("Sex","Age","Osstatus","DFSstatus","OS","DFS","LeftOrRight","LymStatus","NI","VI","Stage","remove","riskscore","riskscoreStatus","ageStatus"))]
#'
#' # cox
#' y <- df[,which(colnames(df) %in% c("Osstatus","OS"))]
#' y <- y[,c(2,1)]
#' colnames(y) <- c("time","status")
#' m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "cox",parallel=F)
#'
#' # binomial
#' y <- df$riskscoreStatus
#' m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "binomial",parallel=F)
#'
#' # gaussian
#' y <- df$riskscore
#' m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F)
#'
#' # if you don't need the result of permutation
#' m<-Lasso.bag(x,y,bootN=3,imputeN=5,imputeN.max = 7,permut.increase = 1,boot.rep = T,a.family = "gaussian",parallel=F,permutation=FALSE)

#   Test Package:              'Cmd + Shift + T'

#' @export

# library(glmnet)
# library(parallel)
# library(POT)
# source("LessPermutation.R")
Lasso.bag <- function(mat,out.mat,bootN=1000,imputeN=1000,imputeN.max=2000,permut.increase=100,boot.rep=TRUE,a.family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),parallel=F,fit.pareto="mle",permutation=TRUE,n.cores=1,rd.seed=89757,plot.freq="full",plot.out=F, use.gpd=F) {
  # bootN is the size of resample sample
  # mat is independent variable
  # out.mat is dependent variable
  # boot.rep is whether :"sampling with return" or not
  # a.family is what kind of regression method to use, it should match the type of out.mat
  # imputeN is the initial permutation times
  # imputeN.max is the max permutation times. Regardless of whether p has meet the requirement
  # permut.increase is: if the initial imputeN times of permutation doesn't meet the requirement, then we add permut.increase times of permutation to get more random values
  # if a.family == "multinomial", then the number of each class of dependent vars should be more than 1
  set.seed(rd.seed)
  rownames(mat) <- c(1:length(rownames(mat)))
  all.num <- bootN * imputeN
  # rownames(out.mat) <- c(1:length(rownames(out.mat)))
  # a simple input judgement and transformation
  out.mat <- as.matrix(out.mat)
  if (a.family!="cox") {
    colnames(out.mat) <- "out"
  }

  # to keep the same type of each column, default is force them to be "numeric"
  zeromat <- matrix(0,length(rownames(mat)),length(colnames(mat)))
  colnames(zeromat) <- colnames(mat)
  rownames(zeromat) <- rownames(mat)
  for (i in colnames(mat)) {
    zeromat[,i] <- as.numeric(mat[,i])
  }
  mat <- zeromat  # it's necessary to substitute in this way for there will be some mistake in later handling if we don't do this

  if (a.family!="cox") {
    zeromat <- matrix(0,length(out.mat),length(colnames(out.mat)))
    colnames(zeromat) <- colnames(out.mat)
  } else {
    zeromat <- matrix(0,length(rownames(out.mat)),length(colnames(out.mat)))
    colnames(zeromat) <- colnames(out.mat)
    rownames(zeromat) <- rownames(out.mat)
  }

  for (j in colnames(out.mat)) {
    zeromat[,j] <- as.numeric(out.mat[,j])
  }
  out.mat <- zeromat

  # simply judge whether dependent vars has the same size of independent vars
  if (a.family!="cox") {
    if(nrow(mat)!=length(out.mat)){  # TODO: i think that it should be dataframe for that will be more adjustable
      # and it allows us to store different types of data in one mat/df
      warning("incoporate length of matrix and outVarianbles, plz check your input ")
      break
    }
  } else {
    if(nrow(mat)!=length(out.mat[,"time"])){  # TODO: i think that it should be dataframe for that will be more adjustable
      # and it allows us to store different types of data in one mat/df
      warning("incoporate length of matrix and outVarianbles, plz check your input ")
      break
    }
  }


  # lasso bag function in individual iteration
  # boot for one time.
  boot.once<-function(index=NULL){  # TODO: change here to check the input type and take actions to exact types
    print(paste(Sys.time(), "I'm boostraping", sep = "--"))
    if(!is.null(index)){
      new.out.mat<-out.mat[index,]  # re-order the dependent variable
    }else{
      new.out.mat<-out.mat
    }
    if (a.family!="cox") {
      runData <- cbind(mat,out=new.out.mat)
    } else {
      runData <- cbind(mat,new.out.mat)
    }
    marker_candidate <- colnames(mat)
    out.vec<-rep(0,length(marker_candidate))
    names(out.vec)<-marker_candidate
    selecVlist1 <- c()

    #index list for lapply
    index.list.bootonce<-list()
    for(i in 1:bootN){
      sampleindex2 <- sample(1:nrow(runData),1*nrow(runData),rep=boot.rep)  # re-sampling, same size
      index.list.bootonce[[i]]<-sampleindex2
    }

    # apply function
    boot.indiv<-function(sampleIdex){
      effectdata <- runData[sampleIdex,]  # the sample we use in permutation, x and y have been matched
      if (a.family!="cox") {
        glmmod<-glmnet(data.matrix(effectdata[,marker_candidate]), effectdata[,"out"], family = a.family)
        cv.glmmod<-cv.glmnet(as.matrix(effectdata[,marker_candidate]), effectdata[,"out"], family = a.family,nfolds = 4)
        print(paste(Sys.time(), "one Lasso finished...", sep = " ---------**---------"))
      } else {
        glmmod<-glmnet(data.matrix(effectdata[,marker_candidate]), effectdata[,c("time","status")], family = a.family)
        cv.glmmod<-cv.glmnet(data.matrix(effectdata[,marker_candidate]), effectdata[,c("time","status")], family = a.family,nfolds = 4)
        print(paste(Sys.time(), "one Lasso finished...", sep = " ---------**---------"))
      }
      best_lambda <- cv.glmmod$lambda.1se
      result<-coef(glmmod, s = best_lambda)
      return(result@Dimnames[[1]][which(as.numeric(result) != 0)])  # return the coef of each feature
    }

    if (parallel == FALSE){
      selecVlist1 <- mclapply(index.list.bootonce, boot.indiv,mc.cores = n.cores)
    } else {
      selecVlist1 <- lapply(index.list.bootonce, boot.indiv)
    }

    tablecount1 <- table(unlist(selecVlist1))
    out.vec[intersect(names(tablecount1), names(out.vec))] <- tablecount1[intersect(names(tablecount1), names(out.vec))]
    print(paste(Sys.time(), "One boostrap finished ...", sep = "--"))
    return(out.vec)  # the output is "what features have been chosen."
  }
  # get observed value
  cat(paste(date(), "", sep=" -- start calculating observed frequency "), '\n')
  observed.fre<-boot.once()  # observed frequency. The true freq

  if (permutation==TRUE) {
    #construct datalist with permutation
    cat(paste(date(), "", sep=" -- permutate index "), '\n')

    get.permutation <- function(N) {
      # N is how many times to do permutations at this function
      # returns out.df
      index.list<-list()
      for (i in 1:N) {
        temp.index <- sample(1:nrow(mat),nrow(mat),rep=F)
        index.list[[i]]<-temp.index
      }  # index.list is of length:imputeN, each index is the random order of the original matrix

      cat(paste(date(), "", sep=" -- permutating "), '\n')
      # do permutation
      if(!parallel){  # multiprocessing
        permut.list<-lapply(index.list,boot.once)
      }else{
        permut.list<-mclapply(index.list,boot.once,mc.cores = n.cores)
        # permut.list<-lapply(index.list,boot.once)
      }

      return(permut.list)
    }

    get.plist <- function(permut.list,N) {
      # permut.list is from permutation
      # N is the total permutation times

      out.df<-do.call(cbind.data.frame, permut.list)  # output as a data frame
      # features saved in df will be shown else the entities will be 0
      colnames(out.df)<-c(1:N)
      cat(paste(date(), "", sep=" -- getting pvalue "), '\n')

      #get permutation pvalue
      pvalue.list<-c()
      add.more <- F
      for(i in 1:length(observed.fre)){
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
                   print(p.value)
                   print("no data is bigger than threshold, we will use traditional p-value")})

        } else {
          good.fit <- "fit_good_enough"
        }
        #if (!exists("p.value")) {
        p.value<-(length(temp.list[temp.list>observed.fre[i]])+1)/N
        #}
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








      return(list(add.more,pvalue.list))
    }


    # TODO: to add GPD here to reduce permutation times.
    # we must first do imputeN times permutation
    permut.list <- get.permutation(imputeN)
    judgement <- get.plist(permut.list,imputeN)
    add.more <- judgement[[1]]

    # TODO: to add more permutation here
    total.imputeN <- imputeN
    if (add.more) {

      new.permut.list <-c(permut.list)
      while (add.more & total.imputeN<imputeN.max) {
        new.perm <- get.permutation(permut.increase)
        new.permut.list <- c(new.permut.list,new.perm)  # up to then all the permut.list
        new.judgement <- get.plist(new.permut.list,total.imputeN)
        add.more <- new.judgement[[1]]
        total.imputeN <- total.imputeN + permut.increase
        pvalue.list <- new.judgement[[2]]
      }
      total.imputeN <- total.imputeN + permut.increase
    } else {
      pvalue.list <- judgement[[2]]
    }

    #FDR calculation
    cat(paste(date(), "", sep=" -- Pvalue adjusting "), '\n')
    # FDR.list<-p.adjust(pvalue.list,method = "bonferroni")
    FDR.list <- p.adjust(pvalue.list,method = "fdr")
    res.df<-data.frame(variate=names(observed.fre),Frequency=observed.fre,P.value=pvalue.list,P.adjust=FDR.list)
    cat(paste(date(), "", sep=" -- Done"), '\n')
  } else {  # do not permutate and output frequency
    res.df<-data.frame(variate=names(observed.fre),Frequency=observed.fre)
    cat(paste(date(), "", sep=" -- Done"), '\n')
  }

  res.df.need <- res.df[res.df$Frequency!=0,]
  if (plot.freq=="part") {
    plot.df <- res.df.need
  } else if (plot.freq=="full") {
    plot.df <- res.df
  } else if (plot.freq!=FALSE & plot.freq!="part" & plot.freq!="full"){
    plot.df <- res.df
    print("Actually you need to set plot.freq correctly, here we will plot all features.")
  }
  if (plot.freq!=FALSE) {
    if (plot.out!=F) {  # for saving files
      pdf(plot.out)
      print(ggplot(plot.df, aes(reorder(variate, -Frequency), Frequency)) +
              geom_bar(stat = "identity") +
              theme_bw() +
              theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
              xlab(label = NULL))
      dev.off()
    }
    # for plot on the screen
    print(ggplot(plot.df, aes(reorder(variate, -Frequency), Frequency)) +
            geom_bar(stat = "identity") +
            theme_bw() +
            theme(axis.text.x  = element_text(angle=45, vjust = 0.9, hjust = 1)) +
            xlab(label = NULL))
  }
  return(res.df)
}
