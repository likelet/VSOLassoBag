
singleFilter<-function(mat,out.mat,additional.info=NULL,filter_method="spearman",filter_thres_method="permutation fdr",filter_thres_fdr=0.05,filter_rank_cutoff=0.05,rd.seed=10867,
                     parallel=F,n.cores=1,FilterInitPermutT=100,FilterIncPermutSt=100,FilterMAXPermutT=2000,permut_result_report=F,
                     silent=F,filter_result_report=T,report_all=T,subdir=""){

  ## permutation_distance_cutoff / estimate_loss_cutoff is the rank threshold
  ## by default, filter_rank_cutoff=0.05, i.e. top 5%; To select more features, use larger filter_rank_cutoff
  ## filter_thres_fdr determines the confidence level of the selected features are actually on the level of the filter_rank_cutoff
  
  time<-0

  COX_simple<-function(index){
    estimate_set<-numeric(p)
    y<-out.mat[index,]  # Permutate or not
    surv.y<-Surv(y[,1],y[,2])
    for (i in 1:p){
        
      # univariable cox regression
      if (is.null(additional.info)){
        cox<-coxph(surv.y~mat[,i])
      }else{
        cox<-coxph(surv.y~mat[,i]+additional.info)
      }
      estimate_set[i]<-cox$coefficients[1]
    }
    names(estimate_set)<-NULL
    return(estimate_set)
  }

  COX_test<-function(index){
    p_set<-numeric(p)
    estimate_set<-numeric(p)
    y<-out.mat[index,]  # Permutate or not
    marker<-colnames(mat)
    surv.y<-Surv(y[,1],y[,2])
    for (i in 1:p){
      # univariable cox regression
      if (is.null(additional.info)){
        cox<-coxph(surv.y~mat[,i])
      }else{
        cox<-coxph(surv.y~mat[,i]+additional.info)
      }
      estimate_set[i]<-cox$coefficients[1]
      p_set[i]<-summary(cox)$sctest[3]  # Logrank Test P-value
    }
    
    names(estimate_set)<-NULL
    names(p_set)<-NULL
    fdr_set<-p.adjust(p_set,method="fdr")
    ro<-order(x=fdr_set)
    ord<-integer(length(ro))
    ord[ro]<-c(1:length(ro))
    report<-data.frame(rank=ord,feature=marker,estimate=estimate_set,P=p_set,P.adjust=fdr_set) # estimate here is actually the coefficient of the feature in the univariate cox model
    return(report)
  }

  cor_simple<-function(index){  ## only return the estimates to compute faster
    estimate_set<-numeric(p)
    y<-out.mat[index]
    for (i in 1:p){
      x<-mat[,i]
      estimate_set[i]<-cor(x,y,method=filter_method)
    }
    return(estimate_set)
  }

  cor_test<-function(index){
    p_set<-numeric(p)
    estimate_set<-numeric(p)
    y<-out.mat[index]
    for (i in 1:p){
      r<-cor.test(mat[,i],y,method=filter_method,exact=FALSE)
      p_set[i]<-r$p.value
      estimate_set[i]<-r$estimate
    }
    fdr_set<-p.adjust(p_set,method="fdr")
    ro<-order(x=abs(estimate_set),decreasing=T)
    ord<-integer(length(ro))
    ord[ro]<-c(1:length(ro))
    report<-data.frame(rank=ord,feature=colnames(mat),estimate=estimate_set,P=p_set,P.adjust=fdr_set)
    return(report)
  }

  jaccard<-function(a,b){  ## Jaccard coefficient
    if (length(a)==0 | length(b)==0){
      return(0)
    }
    a<-unique(a)
    b<-unique(b)
    c<-intersect(a,b)
    d<-union(a,b)
    return((length(c)/length(d)))
  }
  
  rank_similarity_cal<-function(new_rank_list,rank_list){  ## calculate the similarity between two rank lists
    ## new_rank_list is shorter than / partial of the rank_list
    if (length(rank_list)==0){
      return(0)
    }
    d<-c(1:length(rank_list))
    names(d)<-rank_list
    return(cor(d[intersect(new_rank_list,rank_list)],sort(d[intersect(new_rank_list,rank_list)]),method="kendall"))  ## the returned estimate should be >0 otherwise it means getting reversed results between two cycles
  }

  sts<-Sys.time()

  set.seed(rd.seed)
  if (!parallel){
    n.cores<-1
  }
  if (silent){
    filter_result_report<-F
  }
  
  if (!is.null(additional.info)){
    if (class(additional.info)=="vector"){
      additional.info<-as.matrix(additional.info,ncol=1)
    }
  }
  
  p<-ncol(mat)
  n<-nrow(mat)
  sel<-integer(0)
  
  
  if (filter_method!="cox"){
    out.mat<-as.vector(out.mat)
    Ay<-out.mat
    if (filter_method=="spearman"){
      out.mat<-rank(out.mat)
      for (i in 1:ncol(mat)){
        mat[,i]<-rank(mat[,i])
      }
    }
    report<-cor_test(c(1:n))
    estimate_set<-report$estimate
    ro<-order(x=abs(estimate_set),decreasing=T)
  }else{
    # For cox filter, univariable cox regression is applied
    Ay<-out.mat[,1]
    report<-COX_test(c(1:n))
    ro<-order(report$P.adjust)
  }
  
  fdr_sel<-which(report$P.adjust<filter_thres_fdr)
  
  time<-time+difftime(Sys.time(),sts,units="min")
  
  if (filter_thres_method=="traditional fdr"){
    ## Null hypothesis is estimate==0 and alternative hypothesis!=0
	## when using adjusted P value method, it is assumed that original P values are uniformly distributed, e.g. the feature with p==0.05 is also ranked ~5%
    sel<-fdr_sel
    if (filter_result_report){
      g4<-ggplot(report)+geom_point(aes(x=(rank/nrow(report)),y=P.adjust),size=0.3,colour="blue")+
	    geom_point(aes(x=(rank/nrow(report)),y=P),size=0.15,colour="red")+
        ylim(0,filter_thres_fdr)+xlim(0,min(1,max(report$rank[sel]/nrow(report))))+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Variable Rank")+ylab("Adjusted P Value")
      pdf(paste0("Filter Result",subdir,".pdf"))
      print(g4)
      dev.off()
    }
  }
  
  if (filter_thres_method=="rank"){
    ## FDR is ignored in this mode
    sel<-which(report$rank<=filter_rank_cutoff*p)
    if (filter_result_report){
      g4<-ggplot(report)+geom_point(aes(x=(rank/nrow(report)),y=P.adjust),size=0.3,colour="blue")+
        xlim(0,min(1,filter_rank_cutoff))+ylim(0,max(report$P.adjust[sel]))+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Variable Rank")+ylab("Adjusted P Value")
      pdf(paste0("Filter Result",subdir,".pdf"))
      print(g4)
      dev.off()
    }

  }
  
  # filter_thres_method=="permutation fdr" is Dreprecated
  if (filter_thres_method=="permutation fdr" & FALSE){
    
    report$P<-NULL
    report$P.adjust<-NULL
    if (!silent){
      cat("NOTICE: using permutation test to calculate the confidence level of the estimate statistics with designated estimate loss (coefficient change for Cox model) of each variable",'\n')
      if (FilterIncPermutSt<100 | FilterInitPermutT<100){
        cat("WARNING: too small initial permutation times and/or the too small increase step may cause unstable and undesired results",'\n')
      }
    }
    permutation_time<-0
    distance_matrix<-NULL
    feature_selected<-0
    test_done<-F
    max_rank_stable<-F
    sel<-integer(0)
    old_sel<-integer(0)
    sel_cut<-integer(0)
    full_ranklist<-character(0)
    max_sel_rank<-0
    
    while ((permutation_time<FilterMAXPermutT) & !test_done){
      sts<-Sys.time()
      Increase<-0
      if (permutation_time==0){
        Increase<-FilterInitPermutT
      }else{
        Increase<-FilterIncPermutSt
      }
      permutation_time<-permutation_time+Increase
      temp_permut_estimate_matrix<-matrix(nrow=p,ncol=Increase)

      index.list.permutate<-list()
      for (i in 1:Increase){   ## create permutation index for response variable
        permutateindex <- sample(1:length(out.mat),length(out.mat),replace=F)
        index.list.permutate[[i]]<-permutateindex
      }
      time<-time+difftime(Sys.time(),sts,units="min")
      sts<-Sys.time()
      if (parallel & n.cores>1){    ## get estimates for permutation results
        if (filter_method!="cox"){
          res<-mclapply(index.list.permutate, cor_simple,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)
        }else{
          res<-mclapply(index.list.permutate, COX_simple,mc.cores = n.cores,mc.preschedule=TRUE,mc.cleanup=TRUE)
        }
      }else{
        if (filter_method!="cox"){
          res<-lapply(index.list.permutate, cor_simple)
        }else{
          res<-lapply(index.list.permutate, COX_simple)
        }
      }
      time<-time+difftime(Sys.time(),sts,units="min")*n.cores
      sts<-Sys.time()
      for (i in 1:Increase){
        temp_permut_estimate_matrix[,i]<-as.vector(res[[i]])
      }
      temp_distance_matrix<-matrix(nrow=p,ncol=Increase)
      for (i in 1:nrow(temp_distance_matrix)){
        temp_distance_matrix[i,]<-abs(report$estimate[i]-temp_permut_estimate_matrix[i,])
      }
      distance_matrix<-cbind(distance_matrix,temp_distance_matrix)
      matrix_size<-p*permutation_time
      
      rank_matrix<-matrix(nrow=p,ncol=permutation_time)
      rank_matrix[order(distance_matrix,decreasing=T)]<-c(1:length(distance_matrix))/matrix_size
      
      permut_outstand_time<-integer(p)
      permut_sig<-numeric(p)
      
      for (i in 1:p){
        permut_outstand_time[i]<-length(which(rank_matrix[i,]<=filter_rank_cutoff))
      }
      # seems no too much difference for choosing poisson or binomial distribution if permutation time>=100 & filter_rank_cutoff<=0.05
      if (permutation_time>=100 & filter_rank_cutoff<=0.05){
        permut_sig<-ppois(permut_outstand_time,filter_rank_cutoff*permutation_time,lower.tail = F)
      }else{
        permut_sig<-pbinom(permut_outstand_time,prob=filter_rank_cutoff,size=permutation_time,lower.tail = F)
      }
      permut_sig_adjust<-p.adjust(permut_sig,method="fdr")
      
      report$permut_outstand_time<-permut_outstand_time
      report$P<-permut_sig
      report$P.adjust<-permut_sig_adjust
      temp_report<-report
      temp_report$order<-c(1:nrow(temp_report))  ## column order here is the same order as in the report dataframe, so that it is kept after further sort
      temp_report<-temp_report[order(temp_report$P.adjust,temp_report$rank),]  ## sort the temp_report according to P.adjust and rank of |estimate|
      temp_report_selected<-temp_report[which(temp_report$P.adjust<filter_thres_fdr),]
      sel<-temp_report_selected$order
      new_ranklist<-temp_report_selected$rank
      feature_selected<-length(sel)
      
      new_max_sel_rank<-max(temp_report_selected$rank)
      max_rank_diff<-max_sel_rank/new_max_sel_rank
      if (max_rank_diff>1){
        max_rank_diff<-1/max_rank_diff
      }  ##Jaccard Coefficient for selecting features whose rank not larger than the ones already selected according to P.adjust<filter_thres_fdr
      max_rank_diff<-1-max_rank_diff
      max_rank_similarity<-rank_similarity_cal(temp_report$rank[which(temp_report$rank<=new_max_sel_rank)],full_ranklist)
      max_sel_rank<-max(temp_report_selected$rank)  ## must update AFTER the above operation
      
      diff<-1-jaccard(sel,old_sel)
      rank_similarity<-rank_similarity_cal(new_ranklist,full_ranklist)  ## actually compare new rank list of this round and the full rank list of last round
      
      if (!silent){
          cat(paste("Filter Overall Permutation Times:",permutation_time),'\n')
          cat(paste("Features selected by Filter:",feature_selected),'\n')
          cat(paste0("Features selected Difference between Two Cycles (based on Jaccard Coefficient): ",round(diff*100,2),"%"),'\n')
          cat(paste0("Features rank lists similarity (indicated by Kendall Rank Correlation Tau): ",round(rank_similarity,5)),'\n\n')
      }
      
      if ((diff<=5e-5*Increase) & (rank_similarity>=(1-5e-4*Increase))){   ## Consider the permutation test results is stable if the difference of selected feature No. between two cycles is small enough
        # Kendall Tau standard is more strict than Jaccard Coefficient
        # Jaccard Coefficient Difference <= 0.5% & Kendall Tau>=0.95 (95%) for 100 permutation increase step
        test_done<-T
      }else{
        # sav_length<-as.integer(round(length(new_ranklist)*(1-filter_thres_fdr),0))  ## cut off last filter_thres_fdr% of features to increase stability of the results, and try to pass the threshold again
        # sav_length2<-as.integer(round(length(old_sel)*(1-filter_thres_fdr),0))
        # temp_report_selected<-temp_report_selected[c(1:sav_length),]
        # sel_cut<-temp_report_selected$order
        # new_ranklist<-temp_report_selected$rank

        # diff<-1-jaccard(sel_cut,old_sel[c(1:sav_length2)])
        # rank_similarity<-rank_similarity_cal(new_ranklist,full_ranklist)  ## actually compare new rank list of this round and the full rank list of last round
        # feature_selected<-length(sel_cut)
        
        # cat(paste0("Features selected by Filter with last ",round(filter_thres_fdr*100,2),"% cutoff: ",feature_selected),'\n')
        # cat(paste0("Features selected Difference between Two Cycles (based on Jaccard Coefficient) with last ",round(filter_thres_fdr*100,2),"% cutoff: ",round(diff*100,2),"%"),'\n')
        # cat(paste0("Features rank lists similarity (indicated by Kendall Rank Correlation Tau) with last ",round(filter_thres_fdr*100,2),"% cutoff: ",round(rank_similarity,5)),'\n\n')
        
        # if (max_rank_selection){
          # if (max_rank_diff==0){   ## For max rank selection, difference between two rounds need to be 0 to be "stable"
            # test_done<-T
            # max_rank_stable<-T
            # sel<-temp_report$order[which(temp_report$rank<=max_sel_rank)]
          # }
          # cat(paste0("Features selected by Filter with rank not larger than the ones filtered by P.adjust (i.e. max rank selection method): ",new_max_sel_rank),'\n')
          # cat(paste0("Features selected by max rank selection method exceeded the original selected ones: ",round((new_max_sel_rank/feature_selected-1)*100,2),"%"),'\n')
          # cat(paste0("Features selected Difference between Two Cycles (based on Jaccard Coefficient) by max rank selection method:",round(max_rank_diff*100,2),"%"),'\n')
          # cat(paste0("Features rank lists similarity (indicated by Kendall Rank Correlation Tau) by max rank selection method: ",round(max_rank_similarity,5)),'\n\n')
        # }
        
      }
      
      ## update the selected feature set and the rank list; do not move the codes
      old_sel<-sel  ## update before cutoff
      full_ranklist<-temp_report$rank  ## update the rank list
      ## END
      
      if (max_rank_stable & !silent){
        cat(paste0("NOTICE:The stable result is achieved by selecting features with rank not larger than the ones filtered by P.adjust, i.e. max rank selection method."),'\n')
      }
      
      if (permut_result_report){  ## Default is FALSE, use for debug ONLY
        write.table(distance_matrix,paste0(permutation_time,"_DisMatrix.txt"),sep = "\t",col.names = F,row.names = F)
      }
      
      temp_report_selected$order<-NULL
      if (filter_result_report){
        write.table(temp_report_selected,paste0("filter_report_selected",subdir,".txt"),sep="\t",row.names=F,quote=F)  ## output for inspection; selected features only
      }
      
      ## save some memory
      rm(temp_permut_estimate_matrix,temp_distance_matrix,res,index.list.permutate,temp_report,temp_report_selected)
      gc()
      ## done
      time<-time+difftime(Sys.time(),sts,units="min")
      if (FilterIncPermutSt<=0){  ##Increase<=0, run test with initial times and then exit
        if (!silent){
          cat("WARNING: Increasing step <= 0, test with initial times done and now exiting...",'\n')
        }
        break
      }
    }
    
    if (filter_result_report){
    
      yy<-as.vector(distance_matrix)
      selyy<-which(as.vector(rank_matrix)<=filter_rank_cutoff)
      g2<-ggplot()+geom_histogram(aes(x=yy),bins=50,alpha=0.5)+geom_histogram(aes(x=yy[selyy]),bins=50,fill="blue",alpha=0.5)+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Distance / Estimate Loss (Upper Tail Cutoff Showed In Blue)")+ylab("Frequency")
      pdf(paste0("Overall Observed Estimates to Baseline Distance",subdir,".pdf"))
      print(g2)
      dev.off()
    
      max_sel<-max(report$rank[sel])/p+0.01
    
      # g3<-ggplot()+geom_histogram(aes(x=report$permut_outstand_time),bins=50)+
        # theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        # xlab("Estimate Sig Times for Each Variable")+ylab("Freq")
      # pdf("Overall Estimate Sig Times.pdf")
      # print(g3)
      # dev.off()
    
      g4<-ggplot(report)+geom_point(aes(x=(rank/nrow(report)),y=P.adjust),colour="blue",size=0.3)+
        ylim(0,filter_thres_fdr+0.01)+xlim(max(max_sel-0.05,0),max_sel)+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Variable Rank")+ylab("Adjusted P Value")
      pdf(paste0("Filter Result ZoomIn",subdir,".pdf"))
      print(g4)
      dev.off()
      
      g5<-ggplot(report)+geom_point(aes(x=(rank/nrow(report)),y=P.adjust),colour="blue",size=0.3)+
        theme(axis.title = element_text(size=15),axis.text = element_text(size=13))+
        xlab("Variable Rank")+ylab("Adjusted P Value")
      pdf(paste0("Filter Result",subdir,".pdf"))
      print(g5)
      dev.off()
      
      
    }
  }
  
  if (filter_result_report){
    
    # ## heatmap
	# Dreprecated
    # #column order
    # feature_ord<-report[sel,]
    # feature_ord<-feature_ord$feature[order(feature_ord$estimate)]
    # heatmap_matrix<-mat[,feature_ord]
    # #row order
    
    # row_ord<-hclust(dist(heatmap_matrix))$order
    # Ay<-Ay[row_ord]
    # sample_label<-data.frame(text=rep("·",length(Ay)),x=rep(-2,length(Ay)),y=c(1:length(Ay)))
    # heatmap_matrix<-heatmap_matrix[row_ord,]
    # heatmap_matrix<-scale(heatmap_matrix)
    # heatmap_matrix_melt<-reshape2::melt(heatmap_matrix,varnames=c("samples","features"))
    # g6<-ggplot(heatmap_matrix_melt)+geom_tile(aes(x=features,y=samples,fill=value))+geom_text(aes(x=x,y=y,color=Ay,label=text),data=sample_label,size=1)+
    # scale_fill_gradient2(low="#001966", high="#660000", mid="#ffffcc")+scale_color_gradient(low="blue",high="red",name="dependent variable")+
    # theme(axis.title = element_text(size=15),axis.text = element_blank(),axis.ticks=element_blank())
    # pdf(paste0("Heatmap",subdir,".pdf"))
    # print(g6)
    # dev.off()
    
    report_selected<-report[sel,]
    report_selected<-report_selected[order(report_selected$P.adjust,report_selected$rank),]
    write.table(report_selected,paste0("filter_report_selected",subdir,".txt"),sep="\t",row.names=F,quote=F)  ## output for inspection; selected features only
    if (report_all){  ## if the overall variable No. is too large or set by user, only show the results of selected variables instead
      report$selected<-""
      report$selected[sel]<-"*"
      write.table(report[order(report$P.adjust,report$rank),],paste0("filter_report_all",subdir,".txt"),sep="\t",row.names=F,quote=F)  ## output final report, all features, sorted according to P.adjust and then rank of |estimate|
    }
  }
  
  gc()
  if (!silent){
    cat(paste("In total",length(sel),"features were selected by the filter."),'\n')
  }
  re<-list(sel=sel,rank=report$rank[sel],time=time)
  return(re)
}



filters<-function(fmat,fout.mat,silent=FALSE,additional.info=NULL,filter_method="auto",a.family,
                  filter_thres_method="traditional fdr",filter_thres_fdr=0.05,filter_rank_cutoff=0.05,FilterInitPermutT=100,FilterIncPermutSt=100,FilterMAXPermutT=2000,
                  parallel=FALSE,n.cores=1,filter_result_report=TRUE,report_all=TRUE,rd.seed=89757){
  # A Correlation Filter here to reduce the feature size required to analysis, and thus boost the algorithm; parallel boosting computation allowed
  timef<-0
  sel<-integer(0)
  if (is.vector(fout.mat)){  ## Force to be matrix
    fout.mat<-as.matrix(fout.mat,ncol=1)
    colnames(fout.mat)<-"out"
    rownames(fout.mat)<-rownames(fmat)
  }
  
  if (filter_method=="auto"){
    if (a.family=="cox"){
      filter_method<-"cox"
    }else{
      filter_method<-"spearman"
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
  
  if (a.family=="cox" & filter_method!="cox"){
    sav<-which(fout.mat[,2]==1) # only retain patients with event happened
    fmat<-fmat[sav,]
    fout.mat<-as.matrix(fout.mat[sav,1])
  }
  
  if (filter_method!="cox"){
    for (i in 1:ncol(fout.mat)){
      filterre<-singleFilter(mat=fmat,out.mat=fout.mat[,i],filter_method=filter_method,filter_thres_method=filter_thres_method,filter_thres_fdr=filter_thres_fdr,filter_rank_cutoff=filter_rank_cutoff,
                           parallel=parallel,n.cores=n.cores,FilterInitPermutT=FilterInitPermutT,FilterIncPermutSt=FilterIncPermutSt,FilterMAXPermutT=FilterMAXPermutT,additional.info=additional.info,
                           filter_result_report=filter_result_report,report_all=report_all,silent=silent,rd.seed=rd.seed)
      sel<-union(sel,filterre$sel)
      timef<-timef+filterre$time
    }
  }else{
      filterre<-singleFilter(mat=fmat,out.mat=fout.mat,filter_method=filter_method,filter_thres_method=filter_thres_method,filter_thres_fdr=filter_thres_fdr,filter_rank_cutoff=filter_rank_cutoff,
                           parallel=parallel,n.cores=n.cores,FilterInitPermutT=FilterInitPermutT,FilterIncPermutSt=FilterIncPermutSt,FilterMAXPermutT=FilterMAXPermutT,additional.info=additional.info,
                           filter_result_report=filter_result_report,report_all=report_all,silent=silent,rd.seed=rd.seed)
      sel<-union(sel,filterre$sel)
      timef<-timef+filterre$time
  }
  return(list(sel=sel,time=timef))
}
