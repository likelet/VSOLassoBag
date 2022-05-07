# for calculating p value with less permutation.
 # source('gradient_descent.R')
 # source('pareto_simple_opt_func.R')
#' LessPermutation: to reduce permutation times by fitting generalized pareto distribution of the right tail data
#'
#' This is an internal function utilized by VSOLassoBag.
#'
#' @param X a vector of data recording the permutation values
#' @param x0 observed value
#' @param fitting.method method to fit GPD, default is "mle", alternative "gd"(gradient descend)
#' @param search.step the length of step (this param * length(X)) to find threshold. Default 0.01
#' @param fit.cutoff the cutoff of p value to judge whether it fits well to GPD, default is 0.05
#' @param when.to.fit a cutoff to tell how many sample values are bigger than the target value then we don't need to fit GPD. it is a portion.Default 0.05
#' @return p value of the observed value in the permutation test
#' @export
# x <- rgpd(200, 1, 2, 0.25)
# LessPermutation(x,1,fitting.method='gd')

LessPermutation <- function(X, x0, fitting.method="mle",search.step=0.01,fit.cutoff=0.05, when.to.fit=0.05) {
  # using General Pareto Distribution to fit the exceedances and return an estimated p value
  #
  # Args:
  #  X: A Union of input data, e.g. c(1,2,3,4,5,6).
  #  x0: Is the target number you need to be surveillanced
  #  fitting.method: Is the fitting method of General Pareto Distribution(GPD)
  #  search.step: Is the length of step (this param * length(X)) to find threshold
  #  fit.cutoff: the cutoff of p value to judge whether it fits well to GPD
  #  when.to.fit: a cutoff to tell how many sample values are bigger than the target value then we
  #               don't need to fit GPD. it is a portion.
  #
  # Returns:
  #  The estimated p value
  X <- sort(X)
  itertimes <- 0 # number of steps to find threshold
  idx <- 0
  thres <- 0
  goodFit <- "Not_good_enough"
  judgeNexc <- "Pareto"
  thresInRange <- "NO"

  # Sort out those entries beyond thres
  OverThres <- function(X, thres) {
    X <- sort(X)
    out <- c()
    for (x in X) {
      if (x > thres) {
        out <- c(out,x)
      }
    }
    return(out)
  }

  # To find threshold according to steps
  FindThres <- function(X, itertime, search.step) {
    ordered_X <- sort(X)
    idx <- 1 + floor((itertime * search.step) * length(ordered_X))
    if (idx >= length(ordered_X)) {
      onemore <- ordered_X[idx]
    } else {
      onemore <- ordered_X[idx + 1]
    }
    return((ordered_X[idx] + onemore) / 2) #this is accorded to the paper
  }


  k_hyb <- function(theta,X){
    right <- sum(log(1-theta*X));
    left <- -1/length(X);
    return(left*right);
  }

  sigma_hyb <- function(k_hyb,theta){
    return(k_hyb/theta);
  }

  #To estimate k and s in GPD
  adtestGPD.EstKS <- function(X,thres,fitting.method) {
    if(fitting.method=="mgf"){
      return(fitgpd(X, thres, fitting.method, stat="ADR")$param);
    }else if (fitting.method=="gd"){
      theta <- GradientDescent(X)
      # k <- -1 * k_hyb(theta,X)
      k <- k_hyb(theta,X)
      s <- sigma_hyb(k,theta)
      output<-matrix(c(k,s),nrow = 1,ncol = 2)
      output <- as.data.frame(output)
      colnames(output) <- c("shape","scale")
      return(output)
    }else if (fitting.method=="pso"){
      theta <- PSOptTheta(X)
      # k <- -1 * k_hyb(theta,X)
      k <- k_hyb(theta,X)
      s <- sigma_hyb(k,theta)
      output<-matrix(c(k,s),nrow = 1,ncol = 2)
      output <- as.data.frame(output)
      colnames(output) <- c("shape","scale")
      return(output)
    }else{

      return(fitgpd(X, thres, fitting.method)$param);
    }
  }

  #GPD function
  GPD.kZero <- function(z, scale) {
    return(1 - exp(-z / scale))
  }
  GPD.kOthers <- function(z, scale, k) {
    return(1 - (1 - k * z / scale)**(1 / k))
  }

  #Make a list of F(z) for AD Test
  adtestGPD.MakeZ <- function(X,thres, k, s) {
    union.Z <- OverThres(X,thres)
    # print("thres is")
    # print(thres)
    # print(union.Z)
    # k <- as.numeric(EstKS["shape"]) * -1
    # s <- as.numeric(EstKS["scale"])
    z <- union.Z - thres
    if (s < 1e-18) {
      s <- 1e-18
    }
    if (k == 0) {
      return(GPD.kZero(z, s))
    } else {
      return(GPD.kOthers(z, s, k))
    }
  }

  #calculate A^2, which is the estimator of AD Test
  adtestGPD.Asqr <- function(zlist) {
    n <- length(zlist)
    sum.asqr <- 0
    Record.Union <- c(1:n)
    for (idx in Record.Union) {
      left <- (2 * idx - 1) * log(zlist[idx])
      # right <- (2 * n + 1 - 2 * idx) * log(1 - zlist[idx]) # actually it is the same as the downside one
      right <- (2 * idx - 1) * log(1-zlist[n+1-idx])
      total <- left + right
      sum.asqr <- sum.asqr + total
    }
    # adjust A^2 if sample is too small
    if (n<=5) {
      asqr <- -1 * n - sum.asqr / n
      adj.asqr <- asqr * (1 + 0.75 / n + 2.25 / (n * n))
      return(adj.asqr)
    } else {
      return(-1 * n - sum.asqr / n)
    }
  }

  #calculate p value and return whether it fits GPD well or not
  p.record <- function(Asqr) {
    #ref http://www.statisticshowto.com/anderson-darling-test/
    #or more exactly this paper:
    #https://www.jstor.org/stable/1165059?seq=1#page_scan_tab_contents
    if (Asqr >= 0.6) {
      pval <- exp(1.2937 - 5.709 * Asqr + 0.0186 * (Asqr ** 2));
    } else if (Asqr < 0.6 & Asqr > 0.34) {
      pval <- exp(0.9177 - 4.279 * Asqr - 1.38 * Asqr ** 2);
    } else if (Asqr <= 0.34 & Asqr > 0.2) {
      pval <- 1 - exp(-8.318 + 42.796 * Asqr - 59.938 * Asqr ** 2);
    } else if (Asqr <= 0.2) {
      pval <- 1 - exp(-13.436 + 101.14 * Asqr - 223.73 * Asqr ** 2);
    }

    if (pval > fit.cutoff) {
      return("fit_good_enough");
    } else {
      return("Not_good_enough");
    }
  }

  #To see if p meets our requirement
  adtestGPD.IsReach <- function(X, thres, k, s) {
    zlist <- adtestGPD.MakeZ(X,thres, k, s)
    Asqr <- adtestGPD.Asqr(zlist)
    return(p.record(Asqr))
  }

  #count p value
  Calp <- function(X, x0, thres) {
    X <- sort(X)
    N <- length(X)
    Nexc <- length(OverThres(X, thres))
    estks <- adtestGPD.EstKS(X, thres, fitting.method)
    k <- -1 * as.numeric(estks["shape"])
    s <- as.numeric(estks["scale"])
    z <- x0 - thres
    # print(thres)
    if(k == 0){
      return(Nexc * (1 - GPD.kZero(z, s)) / N) # according to paper: "Fewer permutations, more accurate P-values"
    } else {
      return(Nexc * (GPD.kOthers(z, s, k)) / N) # according to paper: "Fewer permutations, more accurate P-values"
    }
  }

  # main function
  while (idx < length(X) & judgeNexc != "Regular_EDF"
         & goodFit == "Not_good_enough") {
    if (length(OverThres(X, x0)) / length(X) >= when.to.fit) {
      judgeNexc <- "Regular_EDF" # using traditional method
    } else {
      thres <- FindThres(X, itertimes, search.step)
      if (length(OverThres(X, thres)) / length(X) >= when.to.fit) {
        idx <- ceiling((1 - when.to.fit) * length(X))
        itertimes <- floor((idx -1) / length(X) / search.step) # find threshold in step length of e.g. 0.01*length(X)
        thres <- FindThres(X, itertimes, search.step)
      }
      judgeNexc <- "Pareto"
      itertimes <- itertimes + 1
      idx <- ceiling(1 + (itertimes * search.step) * length(X)) #Update the idx and itertimes

      #limit k and s or x0 - threshold to acceptable range
      EstKS <- adtestGPD.EstKS(X, thres, fitting.method)
      k <- as.numeric(EstKS["shape"]) * -1
      s <- as.numeric(EstKS["scale"])
      z <- x0 - thres
      if (k > 0 & z < (s/k) & z > 0) {
        # print("situ1")
        thresInRange <- "YES"
        goodFit <- adtestGPD.IsReach(X, thres, k, s)
        itertimes <- itertimes + 1
        idx <- ceiling(1 + (itertimes * search.step) * length(X))
      } else if (k < 0 & z > 0){
        # print("situ2")
        thresInRange <- "YES"
        goodFit <- adtestGPD.IsReach(X, thres, k, s)
        itertimes <- itertimes + 1
        idx <- ceiling(1 + (itertimes * search.step) * length(X))
      } else if (k == 0) {
        # print("situ3")
        thresInRange <- "YES"
        goodFit <- adtestGPD.IsReach(X, thres, k, s)
        itertimes <- itertimes + 1
        idx <- ceiling(1 + (itertimes * search.step) * length(X))
      } else {
        # print("situ4")
        thresInRange <- "NO"
        itertimes <- itertimes + 1
        idx <- ceiling(1 + (itertimes * search.step) * length(X))
      }
    }
  }
  if (judgeNexc == "Pareto") {
    if (thresInRange == "YES") {
      if (goodFit == "fit_good_enough") {
        out.p <- Calp(X, x0, thres)
        return(out.p)
      } else {
        print("add_more_permutation_or_just_use_the_current_p_value")
        return(length(OverThres(X, x0)) / length(X))
      }
    } else {
      print("add_more_permutation_or_just_use_the_current_p_value")
      return(length(OverThres(X, x0)) / length(X))
    }
  } else {
    out.p <- length(OverThres(X, x0)) / length(X)
    return(out.p)
  }
}

adtest.gpd <- function(X,x0,fitting.method="mle",fit.cutoff=0.05){
  if (length(which(X>x0))!=0 & length(X)>=20) {
    OverThres <- function(X, thres) {
      X <- sort(X)
      out <- c()
      for (x in X) {
        if (x > thres) {
          out <- c(out,x)
        }
      }
      return(out)
    }
    find.thres <- function(X,x0){
      X <- sort(X)
      a <- which(X>=x0)[1]
      b <- which(X<=x0)[length(which(X<=x0))]
      return((X[a]+X[b])/2)
    }
    adtestGPD.EstKS <- function(X,thres,fitting.method) {
      if(fitting.method=="mgf"){
        return(fitgpd(X, thres, fitting.method, stat="ADR")$param);
      }else if (fitting.method=="gd"){
        theta <- GradientDescent(X)
        # k <- -1 * k_hyb(theta,X)
        k <- k_hyb(theta,X)
        s <- sigma_hyb(k,theta)
        output<-matrix(c(k,s),nrow = 1,ncol = 2)
        output <- as.data.frame(output)
        colnames(output) <- c("shape","scale")
        return(output)
      }else if (fitting.method=="pso"){
        theta <- PSOptTheta(X)
        # k <- -1 * k_hyb(theta,X)
        k <- k_hyb(theta,X)
        s <- sigma_hyb(k,theta)
        output<-matrix(c(k,s),nrow = 1,ncol = 2)
        output <- as.data.frame(output)
        colnames(output) <- c("shape","scale")
        return(output)
      }else{

        return(fitgpd(X, thres, fitting.method)$param);
      }
    }

    #GPD function
    GPD.kZero <- function(z, scale) {
      return(1 - exp(-z / scale))
    }
    GPD.kOthers <- function(z, scale, k) {
      return(1 - (1 - k * z / scale)**(1 / k))
    }

    #Make a list of F(z) for AD Test
    adtestGPD.MakeZ <- function(X,thres, k, s) {
      union.Z <- OverThres(X,thres)
      # print("thres is")
      # print(thres)
      # print(union.Z)
      # k <- as.numeric(EstKS["shape"]) * -1
      # s <- as.numeric(EstKS["scale"])
      z <- union.Z - thres
      if (s < 1e-18) {
        s <- 1e-18
      }
      if (k == 0) {
        return(GPD.kZero(z, s))
      } else {
        return(GPD.kOthers(z, s, k))
      }
    }

    #calculate A^2, which is the estimator of AD Test
    adtestGPD.Asqr <- function(zlist) {
      n <- length(zlist)
      sum.asqr <- 0
      Record.Union <- c(1:n)
      for (idx in Record.Union) {
        left <- (2 * idx - 1) * log(zlist[idx])
        # right <- (2 * n + 1 - 2 * idx) * log(1 - zlist[idx]) # actually it is the same as the downside one
        right <- (2 * idx - 1) * log(1-zlist[n+1-idx])
        total <- left + right
        sum.asqr <- sum.asqr + total
      }
      # adjust A^2 if sample is too small
      if (n<=5) {
        asqr <- -1 * n - sum.asqr / n
        adj.asqr <- asqr * (1 + 0.75 / n + 2.25 / (n * n))
        return(adj.asqr)
      } else {
        return(-1 * n - sum.asqr / n)
      }
    }

    #calculate p value and return whether it fits GPD well or not
    p.record <- function(Asqr) {
      #ref http://www.statisticshowto.com/anderson-darling-test/
      #or more exactly this paper:
      #https://www.jstor.org/stable/1165059?seq=1#page_scan_tab_contents
      if (Asqr >= 0.6) {
        pval <- exp(1.2937 - 5.709 * Asqr + 0.0186 * (Asqr ** 2));
      } else if (Asqr < 0.6 & Asqr > 0.34) {
        pval <- exp(0.9177 - 4.279 * Asqr - 1.38 * Asqr ** 2);
      } else if (Asqr <= 0.34 & Asqr > 0.2) {
        pval <- 1 - exp(-8.318 + 42.796 * Asqr - 59.938 * Asqr ** 2);
      } else if (Asqr <= 0.2) {
        pval <- 1 - exp(-13.436 + 101.14 * Asqr - 223.73 * Asqr ** 2);
      }

      if (pval > fit.cutoff) {
        return("fit_good_enough");
      } else {
        return("Not_good_enough");
      }
    }
    adtestGPD.IsReach <- function(X, thres, k, s) {
      zlist <- adtestGPD.MakeZ(X,thres, k, s)
      Asqr <- adtestGPD.Asqr(zlist)
      return(p.record(Asqr))
    }
    thres <- find.thres(X,x0)
    estks <- adtestGPD.EstKS(X, thres, fitting.method)
    k <- -1 * as.numeric(estks["shape"])
    s <- as.numeric(estks["scale"])
    reach <- adtestGPD.IsReach(X,thres,k,s)
    return(reach)
  } else {
    return("Not_good_enough")
  }
}
