#permutation GLM beta
#reference from http://jmlr.org/papers/volume16/wang15a/wang15a.pdf
permGLMbeta<-function(mat,out.mat,imputeN=1000,glmthred.p=0.2){
  out.mat.mt<-as.matrix(out.mat)
  mat.mt<-as.matrix(mat)
  cof<-solve(t(mat.mt)%*%mat.mt)%*%t(mat.mt)

  observed.fre <- cof%*%out.mat.mt
get_beta<-function(index,cof,y.out){

  tempy<-y.out[index,]
  return(cof%*%tempy)
}

#permutation index
index.list<-list()
for (i in 1:imputeN) {
  temp.index=sample(1:nrow(mat),nrow(mat),rep=F)
  index.list[[i]]<-temp.index
}

permut.beta.list<-lapply(index.list, get_beta,cof=cof,y.out=out.mat.mt)

out.df<-do.call(cbind.data.frame, permut.beta.list)
colnames(out.df)=c("")
#get permutation pvalue
pvalue.list<-c()
  for(i in 1:length(observed.fre)){
    temp.list<-as.numeric(out.df[row.names(observed.fre)[i],])
    #two tailed
    pvalue.list<-c(pvalue.list,length(temp.list[temp.list>abs(observed.fre[i])])/imputeN)
  }
filter.variable<-observed.fre[pvalue.list<glmthred.p,]

}
