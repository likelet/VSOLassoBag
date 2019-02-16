#' LASSO-bagging: a lasso based variable selecting CART framework

#' @param mat sample matrix that each column represent a variable and rows represent sample data points
#' @param out.mat vector with the same length as the sample size from `mat`
#' @param a.family a string vector
#' @param universe total number of all strings that vec1 and vec2 comes from
#' @return  a P value
#   Test Package:              'Cmd + Shift + T'

#' @export
library(glmnet)
library(parallel)
library(POT)
source("LessPermutation.R")
Lasso.bag <- function(mat,out.mat,bootN=1000,imputeN=1000,imputeN.max=2000,permut.increase=100,boot.rep=TRUE,a.family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),parallel=F,fit.pareto="mle",permutation=TRUE) {
  # bootN is the size of resample sample
  # mat is independent variable
  # out.mat is dependent variable
  # boot.rep is whether :"sampling with return" or not
  # a.family is what kind of regression method to use, it should match the type of out.mat
  # imputeN is the initial permutation times
  # imputeN.max is the max permutation times. Regardless of whether p has meet the requirement
  # permut.increase is: if the initial imputeN times of permutation doesn't meet the requirement, then we add permut.increase times of permutation to get more random values
  # if a.family == "multinomial", then the number of each class of dependent vars should be more than 1
  
  rownames(mat) <- c(1:length(rownames(mat)))
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
    if(!is.null(index)){
      new.out.mat<-out.mat[index,]  # re-order the dependent variable
    }else{
      new.out.mat<-out.mat
    }
    if (a.family!="cox") {
      runData=cbind(mat,out=new.out.mat)
    } else {
      runData=cbind(mat,new.out.mat)
    }
    marker_candidate=colnames(mat)
    out.vec<-rep(0,length(marker_candidate))
    names(out.vec)<-marker_candidate
    selecVlist1=c()
    
    #index list for lapply
    index.list.bootonce<-list()
    for(i in 1:bootN){
      sampleindex2=sample(1:nrow(runData),1*nrow(runData),rep=boot.rep)  # re-sampling, same size
      index.list.bootonce[[i]]<-sampleindex2
    }
    
    # apply function
    boot.indiv<-function(sampleIdex){
      effectdata=runData[sampleIdex,]  # the sample we use in permutation, x and y have been matched
      if (a.family!="cox") {
        glmmod<-glmnet(as.matrix(effectdata[,marker_candidate]), effectdata[,"out"], family = a.family)
        cv.glmmod<-cv.glmnet(as.matrix(effectdata[,marker_candidate]), effectdata[,"out"], family = a.family)
      } else {
        glmmod<-glmnet(as.matrix(effectdata[,marker_candidate]), effectdata[,c("time","status")], family = a.family)
        cv.glmmod<-cv.glmnet(as.matrix(effectdata[,marker_candidate]), effectdata[,c("time","status")], family = a.family)
      }
      best_lambda <- cv.glmmod$lambda.1se
      result<-coef(glmmod, s = best_lambda)
      return(result@Dimnames[[1]][which(result != 0)])  # return the coef of each feature
    }
    
    selecVlist1=lapply(index.list.bootonce, boot.indiv)
    
    tablecount1=table(unlist(selecVlist1))
    out.vec[intersect(names(tablecount1), names(out.vec))] <- tablecount1[intersect(names(tablecount1), names(out.vec))]
    return(out.vec)  # the output is "what features have been chosen." 
  }
  # get observed value
  cat(paste(date(), "", sep=" -- start calculating observed frequency "), '\n')
  observed.fre<-boot.once()  # observed frequency. The true freq
  
  #construct datalist with permutation
  cat(paste(date(), "", sep=" -- permutate index "), '\n')
  index.list<-list()
  
  if (permutation==TRUE) {
    get.permutation <- function(N) {
      # N is how many times to do permutations at this function
      # returns out.df
      for (i in 1:N) {
        temp.index=sample(1:nrow(mat),nrow(mat),rep=F)
        index.list[[i]]<-temp.index
      }  # index.list is of length:imputeN, each index is the random order of the original matrix
      
      cat(paste(date(), "", sep=" -- permutating "), '\n')
      # do permutation
      if(!parallel){  # multiprocessing
        permut.list<-lapply(index.list,boot.once)
      }else{
        permut.list<-mclapply(index.list,boot.once)
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
        tryCatch({p.value <- LessPermutation(as.numeric(as.character(temp.list)),observed.value,fitting.method = fit.pareto)},
                 error=function(e){p.value<-(length(temp.list[temp.list>observed.fre[i]])+1)/N;
                 print("no data is bigger than threshold, we will use traditional p-value")})
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
        good.fit <- adtest.gpd(temp.list,observed.fre[i],fitting.method = fit.pareto)
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
      while (add.more & total.imputeN<=imputeN.max) {
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
    return(res.df)
  } else {  # do not permutate and output frequency
    res.df<-data.frame(variate=names(observed.fre),Frequency=observed.fre)
    cat(paste(date(), "", sep=" -- Done"), '\n')
    return(res.df)
  }
}
