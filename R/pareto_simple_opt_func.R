# The input X must be like c(1,2,3,4,5,6)


inputpretreat <- function(X,threshold){
  need <- X[which(X>threshold)];
  need <- sort(need);
  return(need - threshold)
}

g <- function(theta,X,idx,sumpow){
  #theta is a parameter of this function
  #xi is the X in index i
  #sumpow is -n/sum(ln(1-theta*x))
  xi <- X[idx];
  sumfunc <- 1-theta*xi;
  sumfunc1 <- sumfunc^sumpow;
  result <- 1-sumfunc1;
  return(result);
}

lenX <- function(X){
  #use for count the length of X
  dfx <- as.data.frame(X);
  x_length <- length(rownames(dfx));
  return(x_length);
}

sumpowfunc <- function(theta,X){
  #theta is a parameter of this function
  #X the a series or combination of input
  #for calculating the power number of g(theta)
  sumfunc1 <- sumlogfunc(theta,X);
  x_length <- lenX(X);
  return(-x_length/sumfunc1);
}

sumlogfunc <- function(theta,X){
  #for calculating sum(log(1-theta*x))
  sumfunc1 <- sum(log(1-theta*X));
  return(sumfunc1);
}

sumfunc_gprai_p3 <- function(theta,X){
  #for calculation sum(X/(1-theta*X))
  return(sum(X/(1-theta*X)));
}

gprai_p1 <- function(theta,X,idx,sumpow){
  #for -(1-theta*xi)^sumpow
  xi <- X[idx]
  return(-(1-theta*xi) ** sumpow);
}

gprai_p2 <- function(theta,X,idx,sumlog){
  #sumlog is from sumlogfunc. It is sum(log(1-theta*X))
  xi <- X[idx];
  x_length <- lenX(X)
  a <- x_length*xi;
  b <- 1-theta*xi;
  return(a/(b*sumlog));
}

gprai_p3 <- function(theta,X,idx,sumlog,sum_gprai_p3){
  #return the part 3 of pareto opt function
  #for reducing calculating step, we input this result already calculated.
  xi <- X[idx];
  x_length <- lenX(X);
  a <- log(1-theta*xi);
  b <- sum_gprai_p3;
  c <- sumlog;
  return(x_length*a*b/(c**2));
}

gprai <- function(p1,p2,p3){
  return(p1*(p2-p3));
}

Gpraisum_p <- function(idx,theta,X,sumlog,sumpow,sum_gprai_p3){
  #is a part of func G. for further sum;
  n <- lenX(X);
  gip1 <- gprai_p1(theta = theta,
                   X = X,
                   idx = idx,
                   sumpow = sumpow);
  gip2 <- gprai_p2(theta = theta,
                   X = X,
                   idx = idx,
                   sumlog = sumlog);
  gip3 <- gprai_p3(theta = theta,
                   X=X,
                   idx = idx,
                   sumlog = sumlog,
                   sum_gprai_p3 = sum_gprai_p3);
  giprai <- gprai(p1 = gip1,
              p2 = gip2,
              p3 = gip3);
  gi <- g(theta = theta,
          X = X,
          idx = idx,
          sumpow = sumpow);
  one <- (2*idx-1)*giprai/gi;
  two <- (2*n+1-2*idx)*giprai/(1-gi); 
  
  return(one-two);
}

Gsum_p <- function(idx,theta,X,sumpow){
  #is a part of target function, for sum.
  n <- lenX(X);
  gi <- g(theta = theta,
          X = X,
          idx = idx,
          sumpow = sumpow);
  one <- (2*idx-1)*log(gi);
  two <- (2*n+1-2*idx)*log(1-gi);
  return(one + two)
}

Gprai <- function(theta,X){
  #is the target function's derivatives;
  n <- lenX(X);
  idx_union <- c(1:lenX(X));
  sumlog <- sumlogfunc(theta,X);
  sumpow <- sumpowfunc(theta,X);
  sum_gprai_p3 <- sumfunc_gprai_p3(theta,X);
  
  sum_Gp <- 0;
  for(idx in idx_union){
    givalue <- Gpraisum_p(idx = idx,
                      theta = theta,
                      X = X,
                      sumlog = sumlog,
                      sumpow = sumpow,
                      sum_gprai_p3 = sum_gprai_p3);
    sum_Gp <- sum_Gp + givalue;
  }
  return(-sum_Gp/n);
}

G <- function(theta,X){
  #is the target function to be minimized;
  n <- lenX(X);
  idx_union <- c(1:lenX(X));
  sumpow <- sumpowfunc(theta,X);
  sum_Gp <- 0
  for(idx in idx_union){
    givalue <- Gsum_p(idx = idx,
                      theta = theta,
                      X = X,
                      sumpow = sumpow);
    sum_Gp <- sum_Gp + givalue;
  }
  return(-n-sum_Gp/n);
}


