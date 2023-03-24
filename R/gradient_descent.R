#this script is for gradient_descent
 # source('pareto_simple_opt_func.R')
utils::globalVariables("runif")
derivative <- function(theta,X){
  #X should be like c(1,2,3,4,5)
  return(Gprai(theta = theta,X = X))
}
target_func <- function(theta,X){
  #X should be like c(1,2,3,4,5)
  return(G(theta = theta,X = X))
}


# try to find parameters of function: y = ax^2+bx+c, assume that,
# want to find a and b to minimize sum(axi^2)
# example
# m <- rgpd(200,loc = 0,scale = 1,shape = 0.05)
# The following function is for minimizing the target function of finding parameters for GPD by Gradient Descent
GradientDescent <- function(X,step=0.0001,acc=10e-06,maxIter=500) {
  #step is learning rate
  #acc means acception
  message(X)
  curX <- runif(1,max = 1/max(X))
  foreX <- curX
  iterCount <- 1
  previousGradient <- 1
  while (TRUE) {
    lr <- step / (sqrt(previousGradient + acc))
    gradient <- derivative(curX,X)
    curX <- curX - lr*gradient
    diff <- abs(curX - foreX)
    foreX <- curX
    funcVal <- target_func(curX,X)
    message(paste(c("Optimized at Func = ", funcVal), sep = "", collapse = ""))
    # iterCount <- iterCount + 1
    if (iterCount >= maxIter) {
      optimizedX <- curX
      optimizedVal <- funcVal
      message(paste(c("Optimization converge at Func = ", funcVal),sep = "", collapse = ""))
      break
    }
    if (diff <= acc) {
      optimizedX <- curX
      optimizedVal <- funcVal
      message(paste(c("Optimization converge at Func = ", funcVal),sep = "", collapse = ""))
      # break
      iterCount <- iterCount + 1
    }
  }
  return(optimizedX)
}




