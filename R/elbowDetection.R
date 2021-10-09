## Ref: http://www.icsi.berkeley.edu/pubs/networking/findingakneedle10.pdf

#' Kneedle Algorithm: to detect elbow point(s) on the curve
#'
#' This is an internal function utilized by Lasso.bag.
#'
#' @param res a dataframe with variables and observed frequency
#' @param S numeric, determines how aggressive the elbow points on the curve to be called, smaller means more aggressive and larger means more conservative
#' @param auto.loose if TRUE, will reduce `kneedle.S` in case no elbow point is found with the set `kneedle.S`
#' @param min.S a numeric value determines the minimal value that `kneedle.S` will be loosed to.
#' @param loosing.factor a numeric value range in (0,1), which `kneedle.S` is multiplied by to reduce itself.
#' @return the original input dataframe along with the elbow point indicator "elbow.point" with elbow point(s) marked with "*", "Diff" the difference curve, "Thres" the threshold.
#' @references \href{https://ieeexplore.ieee.org/document/5961514}{Original Kneedle Algorithm}, the algorithm utilized in LassoBag has been modified.
#' @export

kneedle<-function(res,S=10,auto.loose=TRUE,min.S=0.1,loosing.factor=0.5){
  # smoothed spline fitting is applied
  # Input: a data.frame contains data points (x,y), with x and y column specified by "col"; a sensitivity parameter S (smaller means more aggressive and more likely to detect an elbow point, or larger means more conservative)
  # Output: return a vector containing the x of elbow point(s)
  # Normalize, Calculate the Difference, and then pick Elbow point(s)
  MM_normalize<-function(x){
    return((x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))
  }
  reverse_MM_normalize<-function(x,mini,maxi){
    return(x*(maxi-mini)+mini)
  }
  ord<-order(res$Frequency,decreasing=T)
  d<-data.frame(x=integer(nrow(res)),y=res$Frequency)
  d<-d[ord,]
  res<-res[ord,]
  d$x<-c(1:nrow(d))/nrow(d)
  d$y<-max(d$y)-d$y
  miny<-min(d$y,na.rm=TRUE)
  maxy<-max(d$y,na.rm=TRUE)
  d$x<-MM_normalize(d$x)
  d$y<-MM_normalize(d$y)
  D<-data.frame(x=numeric(nrow(d)),y=numeric(nrow(d)))
  for (i in 2:nrow(D)){
    D$x[i]<-d$x[i]-d$x[i-1]
    D$y[i]<-d$y[i]-d$y[i-1]
  }
  candidate<-integer(nrow(D))
  localmin<-integer(nrow(D))
  for (i in 2:(nrow(D)-1)){
    if (D$y[i]>D$y[i-1] & D$y[i]>D$y[i+1]){
      candidate[i]<-1
    }
    if (D$y[i]<D$y[i-1] & D$y[i]<D$y[i+1]){
      localmin[i]<-1
    }
  }
  candidate<-which(candidate==1)
  localmin<-which(localmin==1)
  Dlmx<-D[candidate,]
  Dlmx$i<-candidate
  avediffx<-(d$x[nrow(d)]-d$x[1])/(nrow(d)-1)
  reached<-FALSE
  while (!reached){
    Dlmx$T<-Dlmx$y-S*avediffx
    Dlmx$knee<-FALSE
    
    for (i in 1:nrow(Dlmx)){
      con<-FALSE
      if (i==nrow(Dlmx)){
        r<-nrow(d)
      }else{
        r<-Dlmx$i[i+1]-1  # just before next local maximum
      }
      lm<-localmin[which(localmin %in% c((Dlmx$i[i]+1):r))]
      if (length(lm)>1){
        lm<-lm[1]
      }
      for (j in c((Dlmx$i[i]+1):r)){
        if (D$y[j]<Dlmx$T[i]){
          con<-TRUE
          break
        }
        if (length(lm)>0){
          if (j>=lm){
            break
          }
        }
      }
      Dlmx$knee[i]<-con
    }
    
    knee_x<-Dlmx$i[which(Dlmx$knee)]
    if (length(knee_x)>0){
      reached<-TRUE
      cat(paste0("Using S = ",S," for elbow point dection."),'\n')
    }else{
      if (!auto.loose){
        reached<-TRUE
      }else{
        if (S>=round(S*loosing.factor,5)){
          reached<-TRUE
        }
        S<-round(S*loosing.factor,5)
        if (S<min.S){
          reached<-TRUE
        }
      }
    }
  }
  res$elbow.point<-""
  if (length(knee_x)>0){
    res$elbow.point[knee_x]<-"*"
  }
  res$Diff<-reverse_MM_normalize(D$y,mini=miny,maxi=maxy)
  res$Thres<-reverse_MM_normalize(D$y-S*avediffx,mini=miny,maxi=maxy)
  Th<-0
  for (i in 1:nrow(res)){
    if (res$elbow.point[i]=="*"){
      Th<-res$Thres[i]
    }
    res$Thres[i]<-Th
  }
  return(res)
}

