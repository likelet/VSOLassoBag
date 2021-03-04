


forwardSelection<-function(mat, out.mat, res.df, a.family, additional.info=NULL,local.min.BIC=FALSE){
  ## A step-wise forward selection method
  ## Rank all features (with baggingSelectedFreq>0) according to baggingSelectedFreq, then apply the forward selection method
  ## select the model with optimal (lowest) BIC
  ## when local.min.BIC = TRUE (default), not BIC of all models are calculated, but the BIC of models before reaching the local minimum BIC, and equal number of models after reaching the local minimum, are evaluated; the local minimum BIC is updated through the whole process
  
  ## Input: mat (X matrix), out.mat (Y matrix), res.df (a data.frame)
  ## Output: return res.df along with the BIC info
  
  # BIC_cus<-function(k,n,mse){
    # if (a.family!="cox"){
      # bic<-n*log(mse)+k*log(n)
    # }else{
      # # input log partial likelihood as mse
      # bic<--2*mse+k*log(n)
    # }
    # return(bic)
  # }
  
  cox_bic<-function(mat,out.mat,additional.info=NULL){
    if (!is.matrix(mat)){
      mat<-as.matrix(mat)
    }
    if (is.null(additional.info)){
      x<-mat
      if (ncol(x)<=1){
        return(NA)
      }
      fit<-glmnet(x,out.mat,family="cox",lambda=0)
    }else{
      x<-cbind(additional.info,mat)
      fit<-glmnet(x,out.mat,family="cox",lambda=0)
    }
    bic<-assess.glmnet(fit,x,out.mat,family = "cox")$deviance+ncol(mat)*log(nrow(mat))
    return(bic)
  }
  
  if (is.vector(out.mat)){
    # force out.mat to be a matrix
    out.mat<-as.matrix(out.mat,ncol=1)
    colnames(out.mat)<-"out"
  }
  
  # Rank features first
  ord<-order(res.df$Frequency,decreasing=TRUE)
  res.df<-res.df[ord,]
  res.df$BIC<-NA
  mat<-mat[,ord]
  
  features<-res.df$variate
  n<-nrow(mat)
  
  local.min<-Inf
  early.break<-ncol(mat)
  
  for (k in 1:ncol(mat)){
    print(k)
    x<-as.matrix(mat[,c(1:k)])
    if (a.family!="cox"){
      fit<-lm(out.mat~x)
      bic<-BIC(fit)
      res.df$BIC[k]<-bic
        # # Utilize glmnet function with no penalization to simplify the coding
        # fit<-glmnet(x=mat[,c(1:k)],y=out.mat,family=a.family,lambda=0)
        # mse<-assess.glmnet(fit,newx=mat[,c(1:k)],newy = out.mat,family=a.family)$mse
        # bic<-BIC_cus(k=k,n=n,mse=mse)
        # res.df$BIC[k]<-bic
    }else{
      res.bic<-NA
      res.bic<-cox_bic(x,out.mat,additional.info)
      res.df$BIC[k]<-res.bic
    }
    if (local.min>res.df$BIC[k]){
      local.min<-res.df$BIC[k]
      early.break<-k*2
    }
    if (local.min.BIC & (k>=early.break)){
      break
    }
  }
  res.df$BIC[which(duplicated(res.df$BIC))]<-NA  # ensure using the model with less features when BIC are equal
  return(res.df)
}
