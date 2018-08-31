#' LASSO-bagging: a lasso based variable selecting CART framework

#' @param mat sample matrix that each column represent a variable and rows represent sample data points
#' @param out.mat vector with the same length as the sample size from `mat`
#' @param a.family a string vector
#' @param universe total number of all strings that vec1 and vec2 comes from
#' @return  a P value
#   Test Package:              'Cmd + Shift + T'
#' @export

Lasso.bag <- function(mat,out.mat,bootN=100,imputeN=100,boot.rep=TRUE,a.family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),parallel=F) {
  if(nrow(mat)!=length(out.mat)){
      warning("incoporate length of matrix and outVarianbles, plz check your input ")
      break
  }

  # lasso bag function in individual iteration

  boot.once<-function(index=NULL){
    if(!is.null(index)){
      new.out.mat<-out.mat[index]
    }else{
      new.out.mat<-out.mat
    }
    runData=cbind(mat,out=new.out.mat)
    marker_candidate=colnames(mat)
    out.vec<-rep(0,length(marker_candidate))
    names(out.vec)<-marker_candidate
    selecVlist1=c()
    for(i in 1:bootN){
      sampleindex2=sample(1:nrow(runData),1*nrow(runData),rep=boot.rep)
      effectdata=runData[sampleindex2,]

      glmmod<-glmnet(as.matrix(effectdata[,marker_candidate]), effectdata$out, family = a.family)
      cv.glmmod<-cv.glmnet(as.matrix(effectdata[,marker_candidate]), effectdata$out, family = a.family)
      best_lambda <- cv.glmmod$lambda.1se
      result<-coef(glmmod, s = best_lambda)
      selecVlist1=c(selecVlist1,result@Dimnames[[1]][which(result != 0)])
    }
    tablecount1=table(selecVlist1)
    out.vec[intersect(names(tablecount1), names(out.vec))] <- tablecount1[intersect(names(tablecount1), names(out.vec))]
    return(out.vec)
  }
  # get observed value
  cat(paste(date(), "", sep=" -- start calculating observed frequency "), '\n')
  observed.fre<-boot.once()

  #construct datalist with permutation
  cat(paste(date(), "", sep=" -- permutate index "), '\n')
  index.list<-list()
  for (i in 1:imputeN) {
    temp.index=sample(1:nrow(mat),nrow(mat),rep=F)
    index.list[[i]]<-temp.index
  }

  cat(paste(date(), "", sep=" -- permutating "), '\n')
  # do permutation
  if(!parallel){
    permut.list<-lapply(index.list,boot.once)
  }else{
    permut.list<-mclapply(index.list,boot.once)
  }
  out.df<-do.call(cbind.data.frame, permut.list)

  colnames(out.df)<-c(1:imputeN)
  cat(paste(date(), "", sep=" -- getting pvalue "), '\n')
  #get permutation pvalue
  pvalue.list<-c()
  for(i in 1:length(observed.fre)){
    temp.list<-out.df[names(observed.fre)[i],]
    pvalue.list<-c(pvalue.list,length(temp.list[temp.list>observed.fre[i]])/imputeN)
  }
  #FDR calculation
  cat(paste(date(), "", sep=" -- Pvalue adjusting "), '\n')
  FDR.list<-p.adjust(pvalue.list,method = "bonferroni")
  res.df<-data.frame(variate=names(observed.fre),Frequency=observed.fre,P.value=pvalue.list,P.adjust=FDR.list)
  cat(paste(date(), "", sep=" -- Done"), '\n')
  return(res.df)
}
