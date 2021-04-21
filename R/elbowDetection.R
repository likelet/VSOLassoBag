## Ref: http://www.icsi.berkeley.edu/pubs/networking/findingakneedle10.pdf
# Adapted Needle Algorithm

kneedle<-function(res,S=1){
  # smoothed spline fitting is applied
  # Input: a data.frame contains data points (x,y), with x and y column specified by "col"; a sensitivity parameter S (smaller means more aggressive, or larger means more conservative)
  # Output: return a vector containing the x of elbow point(s)
  # Normalize, Smooth Fitting, Calculate the Difference, and then pick Elbow point(s)
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
  d$x<-c(1:nrow(d))
  d$x<-MM_normalize(d$x)
  d$y<-1-MM_normalize(d$y)
  z<-smooth.spline(d$x,d$y)
  D<-data.frame(x=d$x,y=z$y-d$x)
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
  avediffx<-(d$x[nrow(d)]-d$x[1])/(nrow(d)-1)
  Dlmx$T<-Dlmx$y-S*avediffx
  Dlmx$knee<-F
  
  for (i in 1:nrow(Dlmx)){
    con<-F
    if (i==nrow(Dlmx)){
      for (j in (candidate[i]+1):nrow(D)){
        if (length(intersect(j,localmin))>=1){
        }
        break
        if (D$y[j]<Dlmx$T[i]){
          con<-T
          break
        }
      }
    }else{
      for (j in (candidate[i]+1):(candidate[i+1]-1)){
        if (length(intersect(j,localmin))>=1){
          break
        }
        if (D$y[j]<Dlmx$T[i]){
          con<-T
          break
        }
      }
    }
    Dlmx$knee[i]<-con
  }
  
  knee_x<-Dlmx$x[which(Dlmx$knee)]
  knee_x<-as.integer(round(knee_x*(nrow(d)-1)+1))
  res$fitted_value<-reverse_MM_normalize(1-z$y,mini=min(res$Frequency),maxi=max(res$Frequency))
  res$elbow_point<-""
  if (length(knee_x)>0){
    res$elbow_point[knee_x]<-"*"
  }
  return(res)
}

