utils::globalVariables(c("P.adjust", "P", "pdf", "dev.off", "write.table",
                         "cor.test"))
singleFilter<-function(mat,out.mat,additional.covariable=NULL,
                       filter.method="spearman",filter.thres.method="fdr",filter.thres.P=0.05,filter.rank.cutoff=0.05,
                       filter.min.variables=-Inf,filter.max.variables=Inf,
                       parallel=FALSE,n.cores=1,
                       silent=FALSE,filter.result.report=TRUE,filter.report.all.variables=TRUE,subdir=""){

  ## permutation_distance_cutoff / estimate_loss_cutoff is the rank threshold
  ## by default, filter.rank.cutoff=0.05, i.e. top 5%; To select more features, use larger filter.rank.cutoff
  ## filter.thres.P determines the confidence level of the selected features are actually on the level of the filter.rank.cutoff
  time<-0
  sts<-Sys.time()

  COX_test<-function(index){
    test<-function(index,surv.y,additional.covariable=NULL){
      if (is.null(additional.covariable)){
        cox<-coxph(surv.y~mat[,index])
      }else{
        cox<-coxph(surv.y~mat[,index]+additional.covariable)
      }
      estimate<-cox$coefficients[1]
      p<-summary(cox)$logtest[3]  # Likelihood ratio Test P-value
      names(estimate)<-NULL
      names(p)<-NULL
      return(c(estimate,p))
    }

    p_set<-numeric(p)
    estimate_set<-numeric(p)
    y<-out.mat[index,]  # Permutate or not
    marker<-colnames(mat)
    surv.y<-Surv(y[,1],y[,2])

    if (parallel & n.cores>1){
      gc()
      res <- mclapply(c(1:ncol(mat)), test,surv.y=surv.y,additional.covariable=additional.covariable,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)
    }else{
      res <- lapply(c(1:ncol(mat)), test,surv.y=surv.y,additional.covariable=additional.covariable)
    }
    res<-as.data.frame(res)
    estimate_set<-res[1,]
    p_set<-res[2,]
    # for (i in 1:p){
      # # univariable cox regression
      # if (is.null(additional.covariable)){
        # cox<-coxph(surv.y~mat[,i])
      # }else{
        # cox<-coxph(surv.y~mat[,i]+additional.covariable)
      # }
      # estimate_set[i]<-cox$coefficients[1]
      # p_set[i]<-summary(cox)$logtest[3]  # Likelihood ratio Test P-value
    # }

    names(estimate_set)<-NULL
    names(p_set)<-NULL
    estimate_set<-as.numeric(estimate_set)
    p_set<-as.numeric(p_set)
    fdr_set<-p.adjust(p_set,method="fdr")
    ro<-order(x=fdr_set)
    ord<-integer(length(ro))
    ord[ro]<-c(1:length(ro))
    report<-data.frame(rank=ord,feature=marker,estimate=estimate_set,P=p_set,P.adjust=fdr_set) # estimate here is actually the coefficient of the feature in the univariable cox model
    return(report)
  }

  cor_test<-function(index){
    test<-function(index,y,method){
      res<-cor.test(mat[,index],y,method=method,exact=FALSE)
      estimate<-res$estimate
      names(estimate)<-NULL
      return(c(estimate,res$p.value))
    }
    p_set<-numeric(p)
    estimate_set<-numeric(p)
    y<-out.mat[index]
    if (parallel & n.cores>1){
      gc()
      res <- mclapply(c(1:ncol(mat)), test,y=y,method=filter.method,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)
    }else{
      res <- lapply(c(1:ncol(mat)), test,y=y,method=filter.method)
    }
    res<-as.data.frame(res)
    estimate_set<-res[1,]
    p_set<-res[2,]
    names(estimate_set)<-NULL
    names(p_set)<-NULL
    estimate_set<-as.numeric(estimate_set)
    p_set<-as.numeric(p_set)
    # for (i in 1:p){
      # r<-cor.test(mat[,i],y,method=filter.method,exact=FALSE)
      # p_set[i]<-r$p.value
      # estimate_set[i]<-r$estimate
    # }
    fdr_set<-p.adjust(p_set,method="fdr")
    ro<-order(x=abs(estimate_set),decreasing=TRUE)
    ord<-integer(length(ro))
    ord[ro]<-c(1:length(ro))
    report<-data.frame(rank=ord,feature=colnames(mat),estimate=estimate_set,P=p_set,P.adjust=fdr_set)
    return(report)
  }


 if (silent){
    filter.result.report<-FALSE
  }

  if (!is.null(additional.covariable)){
    if (is(additional.covariable, 'vector')){
      additional.covariable<-as.matrix(additional.covariable,ncol=1)
    }
  }

  p<-ncol(mat)
  n<-nrow(mat)
  sel<-integer(0)


  if (filter.method=="spearman" | filter.method=="pearson" | filter.method=="kendall"){
    out.mat<-as.vector(out.mat)
    if (filter.method=="spearman"){
      out.mat<-rank(out.mat)
      for (i in 1:ncol(mat)){
        mat[,i]<-rank(mat[,i])
      }
    }
    if (filter.method=="kendall"){
      # must be discrete data
      out.mat<-factor(out.mat)
      for (i in 1:ncol(mat)){
        mat[,i]<-factor(mat[,i])
      }
    }
    report<-cor_test(c(1:n))
    estimate_set<-report$estimate
    ro<-order(x=abs(estimate_set),decreasing=TRUE)
  }
  if (filter.method=="cox"){
    # For cox filter, univariable cox regression is applied
    report<-COX_test(c(1:n))
    ro<-order(report$P.adjust,report$P)
  }



  if (filter.thres.method=="fdr"){
    ## Null hypothesis is estimate==0 and alternative hypothesis!=0
    ## when using adjusted P value method, it is assumed that original P values are uniformly distributed, e.g. the feature with p==0.05 is also ranked top ~5%
    sel<-which(report$P.adjust<filter.thres.P)
    if (length(sel)<filter.min.variables){
      sel<-which(report$P.adjust<=sort(report$P.adjust)[filter.min.variables])
    }
    if (length(sel)>filter.max.variables){
      sel<-which(report$P.adjust<=sort(report$P.adjust)[filter.max.variables])
    }
    if (filter.result.report){
      g4<-ggplot(report)+geom_point(aes(x=(rank/nrow(report)),y=P.adjust),size=0.3,colour="blue")+
        geom_point(aes(x=(rank/nrow(report)),y=P),size=0.15,colour="red")+
        ylim(0,filter.thres.P)+xlim(0,min(1,max(report$rank[sel]/nrow(report))))+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Variable Rank")+ylab("Adjusted P Value")
      pdf(paste0("Filter Result",subdir,".pdf"))
      print(g4)
      dev.off()
    }
  }

  if (filter.thres.method=="rank"){
    ## FDR is ignored in this mode
    if (filter.rank.cutoff*p<filter.min.variables){
      sel<-which(report$rank<=filter.min.variables)
    }else{
      if (filter.rank.cutoff*p>filter.max.variables){
        sel<-which(report$rank<=filter.max.variables)
      }else{
        sel<-which(report$rank<=filter.rank.cutoff*p)
      }
    }

    if (filter.result.report){
      g4<-ggplot(report)+geom_point(aes(x=(rank/nrow(report)),y=P.adjust),size=0.3,colour="blue")+
        xlim(0,min(1,filter.rank.cutoff))+ylim(0,max(report$P.adjust[sel]))+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Variable Rank")+ylab("Adjusted P Value")
      pdf(paste0("Filter Result",subdir,".pdf"))
      print(g4)
      dev.off()
    }
  }


  if (filter.result.report){
    report_selected<-report[sel,]
    report_selected<-report_selected[order(report_selected$P.adjust,report_selected$rank),]
    write.table(report_selected,paste0("filter_report_selected",subdir,".txt"),sep="\t",row.names=FALSE,quote=FALSE)  ## output for inspection; selected features only
    if (filter.report.all.variables){  ## if the overall variable No. is too large or set by user, only show the results of selected variables instead
      report$selected<-""
      report$selected[sel]<-"*"
      write.table(report[order(report$P.adjust,report$rank),],paste0("filter_filter.report.all.variables",subdir,".txt"),sep="\t",row.names=FALSE,quote=FALSE)  ## output final report, all features, sorted according to P.adjust and then rank of |estimate|
    }
  }

  gc()
  if (!silent){
    message(paste("In total",length(sel),"features were selected by the filter."),'\n')
  }

  time<-time+difftime(Sys.time(),sts,units="min")
  re<-list(sel=sel,rank=report$rank[sel],time=time)

  #re<-list(sel=sel,rank=report$rank[sel])
  return(re)
}


filters<-function(fmat,fout.mat,additional.covariable=NULL,
                  filter.method="auto",a.family,filter.thres.method="fdr",filter.thres.P=0.05,filter.rank.cutoff=0.05,
                  filter.min.variables=-Inf,filter.max.variables=Inf,
                  parallel=FALSE,n.cores=1,
                  silent=FALSE,filter.result.report=TRUE,filter.report.all.variables=TRUE){
  # A Correlation Filter here to reduce the feature size required to analysis, and thus boost the algorithm; parallel boosting computation allowed
  sel<-integer(0)
  timef<-0
   if (!parallel){
    n.cores<-1
  }


  ## The use of filter.min.variables & filter.max.variables is useful especially in Cox filter
  if (filter.min.variables>ncol(fmat)){
    filter.min.variables<-ncol(fmat)
  }
  if (filter.max.variables>ncol(fmat)){
    filter.max.variables<-ncol(fmat)
  }

  if (is.vector(fout.mat)){  ## Force to be matrix
    fout.mat<-as.matrix(fout.mat,ncol=1)
    colnames(fout.mat)<-"out"
    rownames(fout.mat)<-rownames(fmat)
  }

  if (filter.method=="auto"){
    if (a.family=="cox"){
      filter.method<-"cox"
    }else{
      filter.method<-"spearman"
    }
  }
  if (a.family=="multinomial" & ncol(fout.mat)==1){
    ## transform one-column multinomial dependent variable to multiple dummy variables with 0/1 and apply the filter seperately
    cl<-unique(fout.mat)
    temp_out<-as.vector(fout.mat)  # fout.mat is a one-column matrix, and now transformed to a vector
    temp.fout.mat<-matrix(0,nrow=length(temp_out),ncol=length(cl))
    for (i in 1:length(cl)){
      cur_var<-integer(nrow(fout.mat))
      cur_var[which(temp_out==cl[i])]<-1
      temp.fout.mat[,i]<-cur_var
    }
    fout.mat<-temp.fout.mat
  }

  if (a.family=="cox" & filter.method!="cox"){
    sav<-which(fout.mat[,2]==1) # only retain patients with event happened
    fmat<-fmat[sav,]
    fout.mat<-as.matrix(fout.mat[sav,1])
  }

  if (filter.method!="cox"){
    for (i in 1:ncol(fout.mat)){
      filterre<-singleFilter(mat=fmat,out.mat=fout.mat[,i],filter.method=filter.method,filter.thres.method=filter.thres.method,filter.thres.P=filter.thres.P,filter.rank.cutoff=filter.rank.cutoff,
                             filter.min.variables=filter.min.variables,filter.max.variables=filter.max.variables,
                             parallel=parallel,n.cores=n.cores,additional.covariable=additional.covariable,
                             filter.result.report=filter.result.report,filter.report.all.variables=filter.report.all.variables,silent=silent)
      sel<-union(sel,filterre$sel)
      timef<-timef+filterre$time
    }
  }else{
      filterre<-singleFilter(mat=fmat,out.mat=fout.mat,filter.method=filter.method,filter.thres.method=filter.thres.method,filter.thres.P=filter.thres.P,filter.rank.cutoff=filter.rank.cutoff,
                             filter.min.variables=filter.min.variables,filter.max.variables=filter.max.variables,
                             parallel=parallel,n.cores=n.cores,additional.covariable=additional.covariable,
                             filter.result.report=filter.result.report,filter.report.all.variables=filter.report.all.variables,silent=silent)
      sel<-union(sel,filterre$sel)
      timef<-timef+filterre$time
  }
  return(list(sel=sel,time=timef))
}
