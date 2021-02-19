# provide p-value using the statistical method described in RRLASSO (Park H., et al, 2015; doi:10.1371/journal.pone.0141869)

# Input a dataframe res.df of p*2 (p feature rows with 2 columns, i.e. variate and frequency column; other columns if any, would be neglected)
# Input bootN as bagging times
# feature size is automatically determined as nrow(res.df), i.e. the features included in the res.df

# Return a p-value list

simpleEstimation<-function(res.df,bootN){
  pin<-sum(res.df$Frequency)/(bootN*nrow(res.df))
  pvalue.list<-pbinom(q=res.df$Frequency,size=bootN,prob=pin,lower.tail=F)
  return(pvalue.list)
}

